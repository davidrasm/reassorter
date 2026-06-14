#!/usr/bin/env python3
"""
Simulate ancestral recombination graphs (ARGs) and emit local trees or a
tskit TreeSequence. Two recombination models are available via `mode`:

* "hudson" (default): the coalescent with recombination over a continuous
  genome [0, L); lineages are cut at arbitrary breakpoints.
* "reassortment": a segmented-genome model (as in segmented viruses such as
  influenza). The genome is K discrete unit segments and whole segments are
  reassorted between parents; there is no within-segment recombination.

Model
-----
Going backward in time from the present, two kinds of events compete:

* Coalescence (both modes). When k ancestral lineages remain, each of the
  C(k, 2) pairs coalesces at rate 1 / (ploidy * Ne) per generation, so the
  total coalescence rate is C(k, 2) / (ploidy * Ne). Two lineages chosen
  uniformly at random merge into a common ancestor.

* Recombination (mode="hudson"). Each lineage recombines at rate rho per unit
  of ancestral span it carries; a lineage spanning [left, right) recombines at
  rate rho * (right - left). The lineage splits into two parents at a
  breakpoint drawn uniformly within (left, right): one parent inherits the
  material left of the breakpoint, the other the material to the right.

* Reassortment (mode="reassortment"). Each lineage ancestral for >= 2 segments
  reassorts at the constant rate `reassortment_rate`. At a reassortment every
  segment the lineage carries is assigned independently to the first parent
  with probability `reassortment_bias` (default 0.5) and to the second
  otherwise. Because the two parents are exchangeable, `bias` and 1 - `bias`
  give the same distribution; the bias controls how lopsided splits are
  (0.5 = most even / most reassortment, towards 0 or 1 = mostly no-op events).

The waiting time to the next event is Exponential(total rate over all
lineages). The simulation stops once every genomic position has reached its
most recent common ancestor.

Ancestral material
------------------
Each lineage carries a set of segments [left, right) over which it is
ancestral to the sample, together with the node currently ancestral over each.
A coalescence records edges only where two lineages overlap (the events that
appear in a local tree); once a position is shared by a single lineage it has
found its MRCA and is dropped. This is the succinct (tskit) representation: a
recombination/reassortment is not a node, it is simply a child connecting to
two parents over disjoint intervals. In reassortment mode every segment
boundary is at an integer, so a reassortment splits a lineage along those
boundaries.

Conventions
-----------
* Ne is the effective population size (number of *individuals*).
* ploidy = 1 (default) treats the population as haploid (Ne gene copies,
  pairwise rate 1/Ne). ploidy = 2 treats Ne as diploid (2*Ne gene copies,
  pairwise rate 1/(2*Ne)).
* rho is the per-site, per-generation recombination rate (hudson mode).
* L is the genome length; in reassortment mode it equals the number of
  segments K, with unit segments [0,1), [1,2), ..., [K-1, K).
* Branch lengths are in generations.

Output formats
--------------
* Newick (default): one local/marginal tree per line, each prefixed with the
  genomic interval [left, right) over which it applies (adjacent intervals
  sharing the same tree are merged).
* tskit TreeSequence: the node and edge tables are combined into a
  TreeSequence spanning [0, L). Requires the `tskit` package.
"""

import argparse
import random
import sys
from collections import namedtuple


# Explicit ARG network produced by simulate_arg(record_graph=True). Each node is
# a sample, coalescence, or recombination/reassortment event; each edge carries
# the ancestral segments transmitted up that lineage (parent is the older end).
ARGNode = namedtuple("ARGNode", ("id", "time", "type", "meta"))
ARGEdge = namedtuple("ARGEdge", ("parent", "child", "segments"))
ARGGraph = namedtuple("ARGGraph", ("nodes", "edges"))


class OverlapCounter:
    """Piecewise-constant count of how many lineages cover each genome position.

    Used to detect when a genomic interval has reached its MRCA: a coalescence
    that takes an interval's coverage down to 1 means only one lineage carries
    it, so that interval is fully coalesced and need not be propagated further.
    """

    def __init__(self, L, n):
        self.cells = [[0.0, L, n]]  # contiguous [left, right, count] over [0, L)

    def _index_of(self, x):
        for i, (l, r, _) in enumerate(self.cells):
            if l <= x < r:
                return i
        return len(self.cells) - 1

    def _split_at(self, x):
        if x <= self.cells[0][0] or x >= self.cells[-1][1]:
            return
        i = self._index_of(x)
        l, r, c = self.cells[i]
        if x == l:
            return
        self.cells[i] = [l, x, c]
        self.cells.insert(i + 1, [x, r, c])

    def decrement(self, l, r):
        self._split_at(l)
        self._split_at(r)
        for cell in self.cells:
            if cell[0] >= l and cell[1] <= r:
                cell[2] -= 1

    def count_at(self, x):
        return self.cells[self._index_of(x)][2]

    def breakpoints(self):
        return [cell[0] for cell in self.cells[1:]]


def _node_at(segments, x):
    """Return the node id of the segment of `segments` covering x, else None."""
    for l, r, nid in segments:
        if l <= x < r:
            return nid
    return None


def _coalesce_adjacent(segments):
    """Merge neighbouring segments that carry the same node id."""
    merged = []
    for l, r, nid in segments:
        if merged and merged[-1][2] == nid and merged[-1][1] == l:
            pl, _, pnid = merged[-1]
            merged[-1] = (pl, r, pnid)
        else:
            merged.append((l, r, nid))
    return merged


def _merge_lineages(A, B, time, cov, edges, new_node):
    """Coalesce lineages A and B at `time`, recording edges and returning the
    ancestral segments inherited by the new (parent) lineage.

    A node for the common ancestor is created lazily, only if A and B actually
    overlap somewhere (a coalescence of lineages with disjoint ancestral
    material creates no local-tree node, just a lineage carrying both).
    """
    points = set()
    for l, r, _ in A:
        points.add(l)
        points.add(r)
    for l, r, _ in B:
        points.add(l)
        points.add(r)
    for bp in cov.breakpoints():
        points.add(bp)
    points = sorted(points)

    parent = [None]  # boxed so the closure can assign it once

    def parent_id():
        if parent[0] is None:
            parent[0] = new_node(time, False)
        return parent[0]

    segments = []
    for l, r in zip(points, points[1:]):
        mid = (l + r) / 2.0
        a_node = _node_at(A, mid)
        b_node = _node_at(B, mid)
        if a_node is not None and b_node is not None:
            p = parent_id()
            edges.append((l, r, p, a_node))
            edges.append((l, r, p, b_node))
            cov.decrement(l, r)
            if cov.count_at(mid) >= 2:
                segments.append((l, r, p))   # still shared: keep climbing
            # else: coverage is 1, this interval has found its MRCA -> drop
        elif a_node is not None:
            segments.append((l, r, a_node))
        elif b_node is not None:
            segments.append((l, r, b_node))

    return _coalesce_adjacent(segments)


def _split_lineage(segments, bp):
    """Split a lineage's segments at breakpoint bp into (left, right) parts."""
    left, right = [], []
    for l, r, nid in segments:
        if r <= bp:
            left.append((l, r, nid))
        elif l >= bp:
            right.append((l, r, nid))
        else:
            left.append((l, bp, nid))
            right.append((bp, r, nid))
    return left, right


def _segment_count(segments):
    """Number of unit segments a lineage is ancestral for (integer-aligned)."""
    return sum(r - l for l, r, _ in segments)


def _reassort_split(segments, bias, rng):
    """Split a lineage's segments between two parents by per-segment assortment.

    The genome is integer-aligned (unit segments). Each unit segment [i, i+1)
    the lineage carries is assigned to the first parent with probability `bias`
    and to the second otherwise. Returns (first_parent, second_parent) segment
    lists; if one is empty the event was a no-op.
    """
    first, second = [], []
    for l, r, nid in segments:
        i, end = int(round(l)), int(round(r))
        while i < end:
            target = first if rng.random() < bias else second
            target.append((float(i), float(i + 1), nid))
            i += 1
    return _coalesce_adjacent(first), _coalesce_adjacent(second)


def _weighted_choice(weights, rng):
    """Return an index chosen with probability proportional to `weights`."""
    x = rng.random() * sum(weights)
    cumulative = 0.0
    for i, w in enumerate(weights):
        cumulative += w
        if x < cumulative:
            return i
    return len(weights) - 1


def simulate_arg(n, Ne, L, ploidy=1, rng=None, *, mode="hudson", rho=0.0,
                 reassortment_rate=0.0, reassortment_bias=0.5,
                 record_graph=False):
    """Simulate one ARG for `n` samples over a genome of length `L`.

    Parameters
    ----------
    n : int
        Number of sampled individuals (>= 1).
    Ne : float
        Effective population size (number of individuals, > 0).
    L : float
        Genome length (> 0). In reassortment mode this is the number of
        segments K, with unit segments [0,1), ..., [K-1, K).
    ploidy : int
        1 = haploid (default), 2 = diploid. Sets gene copies = ploidy * Ne.
    rng : random.Random, optional
        Random source (pass one with a fixed seed for reproducibility).
    mode : {"hudson", "reassortment"}
        Recombination model (keyword-only).
    rho : float
        [hudson] per-site, per-generation recombination rate (>= 0).
    reassortment_rate : float
        [reassortment] per-lineage reassortment rate (>= 0); applies only to
        lineages ancestral for >= 2 segments.
    reassortment_bias : float
        [reassortment] probability each segment goes to the first parent,
        in [0, 1] (default 0.5).
    record_graph : bool
        If True, also build and return an explicit ARGGraph containing every
        coalescence and recombination/reassortment node (see Returns).

    Returns
    -------
    node_time : list of float
        Time of each node, indexed by node id (samples are ids 0..n-1).
    node_is_sample : list of bool
    edges : list of (left, right, parent, child)
        One record per branch; each spans the interval over which the
        parent-child relationship holds.
    events : list
        One entry per realized recombination/reassortment event (hudson:
        breakpoint positions; reassortment: the index of the lineage split).
    graph : ARGGraph or None
        The explicit network of sample, coalescence, and recombination nodes
        if `record_graph` is True, else None. Each edge carries the ancestral
        segments transmitted up that lineage.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if Ne <= 0:
        raise ValueError("Ne must be > 0")
    if L <= 0:
        raise ValueError("L must be > 0")
    if ploidy < 1:
        raise ValueError("ploidy must be >= 1")
    if mode not in ("hudson", "reassortment"):
        raise ValueError("mode must be 'hudson' or 'reassortment'")
    if mode == "hudson" and rho < 0:
        raise ValueError("rho must be >= 0")
    if mode == "reassortment":
        if reassortment_rate < 0:
            raise ValueError("reassortment_rate must be >= 0")
        if not 0.0 <= reassortment_bias <= 1.0:
            raise ValueError("reassortment_bias must be in [0, 1]")

    rng = rng or random.Random()
    pair_rate = 1.0 / (ploidy * Ne)

    node_time, node_is_sample = [], []

    def new_node(time, is_sample):
        node_time.append(time)
        node_is_sample.append(is_sample)
        return len(node_time) - 1

    # Each lineage is a list of (left, right, node_id) segments. Every sample
    # starts as a lineage carrying the whole genome at the present (t = 0).
    pool = [[(0.0, L, new_node(0.0, True))] for _ in range(n)]

    # Optional explicit ARG network. `lineage_node[i]` is the graph node at the
    # bottom of pool[i]'s current upward stretch, kept in lockstep with `pool`.
    graph_nodes = lineage_node = None
    graph_edges = []
    if record_graph:
        graph_nodes = [ARGNode(i, 0.0, "sample", None) for i in range(n)]
        lineage_node = list(range(n))

    def new_graph_node(time, ntype, meta=None):
        graph_nodes.append(ARGNode(len(graph_nodes), time, ntype, meta))
        return len(graph_nodes) - 1

    cov = OverlapCounter(L, n)
    edges = []
    events = []
    t = 0.0

    while len(pool) > 1:
        k = len(pool)
        coal_rate = k * (k - 1) / 2.0 * pair_rate

        # Per-lineage weights for the non-coalescent event and its total rate.
        if mode == "hudson":
            weights = [seg[-1][1] - seg[0][0] for seg in pool]   # ancestral span
            event_rate = rho * sum(weights)
        else:  # reassortment: only lineages with >= 2 segments can reassort
            weights = [reassortment_rate if _segment_count(seg) >= 2 else 0.0
                       for seg in pool]
            event_rate = sum(weights)

        total_rate = coal_rate + event_rate
        if total_rate <= 0.0:
            break
        t += rng.expovariate(total_rate)

        if rng.random() * total_rate < coal_rate:
            i, j = rng.sample(range(k), 2)
            A, B = pool[i], pool[j]
            merged = _merge_lineages(A, B, t, cov, edges, new_node)
            if record_graph:
                cgid = new_graph_node(t, "coalescence")
                graph_edges.append(ARGEdge(cgid, lineage_node[i], list(A)))
                graph_edges.append(ARGEdge(cgid, lineage_node[j], list(B)))
            for idx in sorted((i, j), reverse=True):
                pool.pop(idx)
                if record_graph:
                    lineage_node.pop(idx)
            if merged:
                pool.append(merged)
                if record_graph:
                    lineage_node.append(cgid)
            continue

        # Otherwise a recombination (hudson) or reassortment event.
        li = _weighted_choice(weights, rng)
        lineage = pool[li]
        if mode == "hudson":
            bp = rng.uniform(lineage[0][0], lineage[-1][1])
            left_part, right_part = _split_lineage(lineage, bp)
            if record_graph:
                rgid = new_graph_node(t, "recombination", {"breakpoint": bp})
                graph_edges.append(
                    ARGEdge(rgid, lineage_node[li], list(lineage)))
            pool.pop(li)
            pool.extend((left_part, right_part))
            if record_graph:
                lineage_node.pop(li)
                lineage_node.extend((rgid, rgid))
            events.append(bp)
        else:
            part_a, part_b = _reassort_split(lineage, reassortment_bias, rng)
            if part_a and part_b:        # both parents inherit material
                if record_graph:
                    meta = {"to_first": [(l, r) for l, r, _ in part_a],
                            "to_second": [(l, r) for l, r, _ in part_b]}
                    rgid = new_graph_node(t, "recombination", meta)
                    graph_edges.append(
                        ARGEdge(rgid, lineage_node[li], list(lineage)))
                pool.pop(li)
                pool.extend((part_a, part_b))
                if record_graph:
                    lineage_node.pop(li)
                    lineage_node.extend((rgid, rgid))
                events.append(li)
            # else: every segment landed on one parent -> a no-op event

    graph = ARGGraph(graph_nodes, graph_edges) if record_graph else None
    return node_time, node_is_sample, edges, events, graph


def squash_edges(edges):
    """Merge edges that share a parent and child across adjacent intervals."""
    edges = sorted(edges, key=lambda e: (e[2], e[3], e[0]))
    out = []
    for l, r, p, c in edges:
        if out and out[-1][2] == p and out[-1][3] == c and out[-1][1] == l:
            pl, _, pp, pc = out[-1]
            out[-1] = (pl, r, pp, pc)
        else:
            out.append((l, r, p, c))
    return out


def _leaf_name(node, node_is_sample):
    return f"n{node + 1}" if node_is_sample[node] else f"i{node}"


def _build_newick(active_edges, node_time, node_is_sample, decimals, label_internal):
    children = {}
    child_nodes = set()
    for p, c in active_edges:
        children.setdefault(p, []).append(c)
        child_nodes.add(c)
    root = next(p for p in children if p not in child_nodes)

    def rec(node):
        if node not in children:
            return _leaf_name(node, node_is_sample)
        parts = []
        for child in children[node]:
            length = node_time[node] - node_time[child]
            parts.append(f"{rec(child)}:{length:.{decimals}f}")
        inner = "(" + ",".join(parts) + ")"
        if label_internal:
            inner += _leaf_name(node, node_is_sample)
        return inner

    return rec(root) + ";"


def marginal_trees(node_time, node_is_sample, edges, L,
                   decimals=6, label_internal=False):
    """Return the local/marginal trees of the ARG as (left, right, newick).

    Adjacent intervals whose tree is identical are merged into one entry.
    """
    points = {0.0, L}
    for l, r, _, _ in edges:
        points.add(l)
        points.add(r)
    points = sorted(p for p in points if 0.0 <= p <= L)

    trees = []
    prev_key = None
    for a, b in zip(points, points[1:]):
        if a >= b:
            continue
        active = [(p, c) for l, r, p, c in edges if l <= a and r >= b]
        key = frozenset(active)
        if key == prev_key:
            la, _, nwk = trees[-1]
            trees[-1] = (la, b, nwk)
            continue
        nwk = _build_newick(active, node_time, node_is_sample,
                            decimals, label_internal)
        trees.append((a, b, nwk))
        prev_key = key
    return trees


def to_tree_sequence(node_time, node_is_sample, edges, L):
    """Build a tskit TreeSequence from the ARG tables. Requires `tskit`."""
    try:
        import tskit
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "tskit is required for TreeSequence output. Install it with "
            "`pip install tskit`."
        ) from exc

    tables = tskit.TableCollection(sequence_length=L)
    for is_sample, time in zip(node_is_sample, node_time):
        flags = tskit.NODE_IS_SAMPLE if is_sample else 0
        tables.nodes.add_row(flags=flags, time=time)
    for left, right, parent, child in squash_edges(edges):
        tables.edges.add_row(left=left, right=right, parent=parent, child=child)
    tables.sort()
    return tables.tree_sequence()


def simplify_arg(graph):
    """Return a copy of `graph` with degree-2 nodes suppressed.

    Any non-sample node with exactly one distinct parent and one distinct child
    carries no branching information; it is removed and its two incident edges
    are spliced into one. This also collapses 'bubbles' where a recombination is
    immediately undone by the two recombinant lineages coalescing back together
    (a recombination node and coalescence node joined by a parallel edge pair).
    Suppression is repeated until no such node remains.
    """
    nodes = {nd.id: nd for nd in graph.nodes}
    alive = set(nodes)

    # Edges keyed by (parent, child); parallel edges (bubbles) merge into one.
    seg = {}
    for edge in graph.edges:
        key = (edge.parent, edge.child)
        seg[key] = _coalesce_adjacent(sorted(seg.get(key, []) + list(edge.segments)))

    def parents_of(v):
        return {p for (p, c) in seg if c == v}

    def children_of(v):
        return {c for (p, c) in seg if p == v}

    changed = True
    while changed:
        changed = False
        for v in list(alive):
            if nodes[v].type == "sample":
                continue
            par, chi = parents_of(v), children_of(v)
            if len(par) == 1 and len(chi) == 1:
                p, c = next(iter(par)), next(iter(chi))
                if p == c:
                    continue
                spliced = seg.get((v, c), [])
                for key in [k for k in seg if v in k]:
                    del seg[key]
                seg[(p, c)] = _coalesce_adjacent(
                    sorted(seg.get((p, c), []) + list(spliced)))
                alive.discard(v)
                changed = True

    new_nodes = [nodes[v] for v in alive]
    new_edges = [ARGEdge(p, c, s) for (p, c), s in seg.items()]
    return ARGGraph(new_nodes, new_edges)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Simulate ancestral recombination graphs (ARGs) under the "
                    "coalescent with recombination or reassortment.")
    p.add_argument("-n", "--num-samples", type=int, required=True,
                   help="number of sampled individuals")
    p.add_argument("-Ne", "--Ne", type=float, required=True,
                   help="effective population size (number of individuals)")
    p.add_argument("--mode", choices=("hudson", "reassortment"),
                   default="hudson",
                   help="recombination model (default: hudson)")
    p.add_argument("--rho", type=float, default=0.0,
                   help="[hudson] per-site, per-generation recombination rate")
    p.add_argument("-L", "--genome-length", type=float, default=None,
                   help="[hudson] genome length")
    p.add_argument("--segments", type=int, default=None,
                   help="[reassortment] number of genome segments K")
    p.add_argument("--reassortment-rate", type=float, default=0.0,
                   help="[reassortment] per-lineage reassortment rate")
    p.add_argument("--reassortment-bias", type=float, default=0.5,
                   help="[reassortment] P(segment -> first parent), default 0.5")
    p.add_argument("--ploidy", type=int, default=1, choices=(1, 2),
                   help="1 = haploid (default), 2 = diploid")
    p.add_argument("--format", choices=("newick", "tskit"), default="newick",
                   help="output format (default: newick)")
    p.add_argument("-r", "--replicates", type=int, default=1,
                   help="number of independent ARGs to simulate (default: 1)")
    p.add_argument("--seed", type=int, default=None,
                   help="random seed for reproducibility")
    p.add_argument("--decimals", type=int, default=6,
                   help="decimal places for branch lengths (default: 6)")
    p.add_argument("--label-internal", action="store_true",
                   help="include internal node labels in the Newick output")
    p.add_argument("-o", "--output", default="-",
                   help="output file ('-' for stdout, the default)")
    p.add_argument("--plot", default=None,
                   help="draw the ARG network to this image file (e.g. arg.png) "
                        "instead of writing trees")
    p.add_argument("--no-simplify", dest="simplify", action="store_false",
                   help="[--plot] keep degree-2 nodes instead of simplifying "
                        "the ARG before drawing")
    p.set_defaults(simplify=True)
    args = p.parse_args(argv)

    if args.mode == "hudson":
        if args.genome_length is None:
            p.error("hudson mode requires -L/--genome-length")
        genome_length = args.genome_length
    else:
        if args.segments is None:
            p.error("reassortment mode requires --segments")
        if args.segments < 1:
            p.error("--segments must be >= 1")
        genome_length = float(args.segments)

    record_graph = args.plot is not None
    rng = random.Random(args.seed)
    sims = [simulate_arg(args.num_samples, args.Ne, genome_length,
                         ploidy=args.ploidy, rng=rng, mode=args.mode,
                         rho=args.rho,
                         reassortment_rate=args.reassortment_rate,
                         reassortment_bias=args.reassortment_bias,
                         record_graph=record_graph)
            for _ in range(args.replicates)]

    if args.plot is not None:
        from viz import draw_arg
        for rep, (_nt, _ns, _edges, _events, graph) in enumerate(sims):
            if args.replicates == 1:
                path = args.plot
            else:
                base, dot, ext = args.plot.rpartition(".")
                path = f"{base}.{rep}.{ext}" if dot else f"{args.plot}.{rep}"
            draw_arg(graph, mode=args.mode, save=path, simplify=args.simplify)
            sys.stdout.write(f"wrote {path}\n")
        return

    event_label = "recombination" if args.mode == "hudson" else "reassortment"

    if args.format == "newick":
        out = sys.stdout if args.output == "-" else open(args.output, "w")
        try:
            for rep, (ntime, nsample, edges, events, _graph) in enumerate(sims):
                trees = marginal_trees(ntime, nsample, edges, genome_length,
                                       decimals=args.decimals,
                                       label_internal=args.label_internal)
                if args.replicates > 1:
                    out.write(f"# replicate {rep}\n")
                out.write(f"# {len(events)} {event_label} event(s), "
                          f"{len(trees)} local tree(s)\n")
                for left, right, nwk in trees:
                    lo = f"{left:.{args.decimals}f}"
                    hi = f"{right:.{args.decimals}f}"
                    out.write(f"[{lo}, {hi}): {nwk}\n")
        finally:
            if out is not sys.stdout:
                out.close()
        return

    # tskit format
    for rep, (ntime, nsample, edges, _events, _graph) in enumerate(sims):
        ts = to_tree_sequence(ntime, nsample, edges, genome_length)
        if args.output == "-":
            sys.stdout.write(
                f"TreeSequence {rep}: {ts.num_samples} samples, "
                f"{ts.num_nodes} nodes, {ts.num_edges} edges, "
                f"{ts.num_trees} tree(s), sequence_length={ts.sequence_length}\n")
        else:
            if args.replicates == 1:
                path = args.output
            else:
                base, dot, ext = args.output.rpartition(".")
                path = f"{base}.{rep}.{ext}" if dot else f"{args.output}.{rep}"
            ts.dump(path)
            sys.stdout.write(f"wrote {path}\n")


if __name__ == "__main__":
    main()
