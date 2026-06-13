#!/usr/bin/env python3
"""
Simulate genealogies (gene trees) under Kingman's coalescent and emit Newick.

Model
-----
Going backward in time from the present, when k ancestral lineages remain,
each of the C(k, 2) pairs coalesces at rate 1 / (number of gene copies) per
generation. The total coalescence rate is therefore

    rate(k) = C(k, 2) / (ploidy * Ne)

and the waiting time to the next coalescence is Exponential(rate(k)),
i.e. mean (ploidy * Ne) / C(k, 2) generations. At each event two lineages,
chosen uniformly at random, merge into a common ancestor.

Conventions
-----------
* Ne is the effective population size (number of *individuals*).
* ploidy = 1 (default) treats the population as haploid (Ne gene copies,
  pairwise rate 1/Ne, E[T_MRCA] ~ 2*Ne*(1 - 1/n) generations).
* ploidy = 2 treats Ne as diploid, so there are 2*Ne gene copies and the
  pairwise rate is 1 / (2*Ne), giving E[T_MRCA] ~ 4*Ne*(1 - 1/n) generations.

Branch lengths in the output are in generations.

Output formats
--------------
* Newick (default): one tree per line.
* tskit TreeSequence: every node (sampled leaves and coalescent ancestors) is
  recorded in a node table with its time, and every branch is recorded in an
  edge table; the two tables are combined into a TreeSequence once the sample's
  MRCA has been reached. Requires the `tskit` package.
"""

import argparse
import random
import sys
from itertools import count


class Node:
    """A node in the genealogy. `time` is height: generations before present."""

    __slots__ = ("name", "time", "children")

    def __init__(self, name=None, time=0.0, children=None):
        self.name = name
        self.time = time
        self.children = children if children is not None else []

    @property
    def is_leaf(self):
        return not self.children


def simulate_coalescent(n, Ne, ploidy=1, rng=None):
    """
    Simulate one coalescent genealogy for `n` samples.

    Parameters
    ----------
    n : int
        Number of sampled individuals (>= 1).
    Ne : float
        Effective population size (number of individuals, > 0).
    ploidy : int
        1 = haploid (default), 2 = diploid. Sets gene copies = ploidy * Ne.
    rng : random.Random, optional
        Random source (pass one with a fixed seed for reproducibility).

    Returns
    -------
    Node
        The root of the genealogy (the MRCA), or the single leaf if n == 1.
    """
    if n < 1:
        raise ValueError("n must be >= 1")
    if Ne <= 0:
        raise ValueError("Ne must be > 0")
    if ploidy < 1:
        raise ValueError("ploidy must be >= 1")

    rng = rng or random.Random()
    gene_copies = ploidy * Ne
    pair_rate = 1.0 / gene_copies

    # Active lineages, each starting as a sampled leaf at the present (t = 0).
    lineages = [Node(name=f"n{i}", time=0.0) for i in range(1, n + 1)]
    t = 0.0
    internal_id = count(1)

    while len(lineages) > 1:
        k = len(lineages)
        total_rate = (k * (k - 1) / 2.0) * pair_rate
        t += rng.expovariate(total_rate)

        # Pick two distinct lineages uniformly at random and merge them.
        i, j = rng.sample(range(k), 2)
        parent = Node(name=f"i{next(internal_id)}", time=t,
                      children=[lineages[i], lineages[j]])
        for idx in sorted((i, j), reverse=True):
            lineages.pop(idx)
        lineages.append(parent)

    return lineages[0]


def to_newick(root, decimals=6, label_internal=False):
    """Render a genealogy rooted at `root` as a Newick string."""

    def fmt(x):
        return f"{x:.{decimals}f}"

    def rec(node):
        if node.is_leaf:
            return node.name
        parts = []
        for child in node.children:
            branch_length = node.time - child.time
            parts.append(f"{rec(child)}:{fmt(branch_length)}")
        inner = "(" + ",".join(parts) + ")"
        if label_internal and node.name:
            inner += node.name
        return inner

    return rec(root) + ";"


def _tree_to_records(root, sequence_length=1.0):
    """
    Flatten a genealogy into tskit-style node and edge records.

    This is the tskit-independent core of the TreeSequence conversion: it walks
    the genealogy once and produces plain-Python records that mirror exactly
    what gets written into the tskit node and edge tables. Keeping it separate
    means the table contents can be validated without tskit installed.

    Returns
    -------
    node_rows : list of (is_sample: bool, time: float)
        Indexed by node id. Sampled leaves come first (ids 0..n-1), ordered by
        their sample number, followed by the coalescent (internal) nodes.
    edge_rows : list of (left, right, parent_id, child_id)
        One row per branch. Every branch spans the whole genome [0, L) because
        without recombination there is a single tree along the sequence.
    """
    leaves, internals = [], []

    def visit(node):
        if node.is_leaf:
            leaves.append(node)
        else:
            for child in node.children:
                visit(child)
            internals.append(node)  # post-order: a node follows its children

    visit(root)

    def sample_number(node):
        # Leaf names are 'n<k>'; order samples so tskit ids 0..n-1 track labels.
        try:
            return int(node.name[1:])
        except (TypeError, ValueError):
            return float("inf")

    leaves.sort(key=sample_number)

    node_id = {}
    node_rows = []
    for leaf in leaves:                       # samples first -> ids 0..n-1
        node_id[leaf] = len(node_rows)
        node_rows.append((True, leaf.time))
    for node in internals:                    # coalescent ancestors next
        node_id[node] = len(node_rows)
        node_rows.append((False, node.time))

    edge_rows = []
    for node in internals:
        parent = node_id[node]
        for child in node.children:
            edge_rows.append((0.0, sequence_length, parent, node_id[child]))

    return node_rows, edge_rows


def to_tree_sequence(root, sequence_length=1.0):
    """
    Convert a genealogy rooted at `root` into a tskit TreeSequence.

    Each node (sampled leaf or coalescent ancestor) is added to the node table
    with its time; each branch is added to the edge table; the tables are then
    sorted and turned into a TreeSequence. Requires the `tskit` package.
    """
    try:
        import tskit
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise ImportError(
            "tskit is required for TreeSequence output. Install it with "
            "`pip install tskit`."
        ) from exc

    node_rows, edge_rows = _tree_to_records(root, sequence_length)

    tables = tskit.TableCollection(sequence_length=sequence_length)
    for is_sample, t in node_rows:
        flags = tskit.NODE_IS_SAMPLE if is_sample else 0
        tables.nodes.add_row(flags=flags, time=t)
    for left, right, parent, child in edge_rows:
        tables.edges.add_row(left=left, right=right, parent=parent, child=child)

    tables.sort()              # tskit requires edges sorted by parent time
    return tables.tree_sequence()


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Simulate coalescent genealogies and output Newick trees.")
    p.add_argument("-n", "--num-samples", type=int, required=True,
                   help="number of sampled individuals")
    p.add_argument("-Ne", "--Ne", type=float, required=True,
                   help="effective population size (number of individuals)")
    p.add_argument("--ploidy", type=int, default=1, choices=(1, 2),
                   help="1 = haploid (default), 2 = diploid")
    p.add_argument("--format", choices=("newick", "tskit"), default="newick",
                   help="output format (default: newick)")
    p.add_argument("--sequence-length", type=float, default=1.0,
                   help="genome length for tskit edges (default: 1.0)")
    p.add_argument("-r", "--replicates", type=int, default=1,
                   help="number of independent trees to simulate (default: 1)")
    p.add_argument("--seed", type=int, default=None,
                   help="random seed for reproducibility")
    p.add_argument("--decimals", type=int, default=6,
                   help="decimal places for branch lengths (default: 6)")
    p.add_argument("--label-internal", action="store_true",
                   help="include internal node labels in the Newick output")
    p.add_argument("-o", "--output", default="-",
                   help="output file ('-' for stdout, the default)")
    args = p.parse_args(argv)

    rng = random.Random(args.seed)
    roots = [simulate_coalescent(args.num_samples, args.Ne,
                                 ploidy=args.ploidy, rng=rng)
             for _ in range(args.replicates)]

    if args.format == "newick":
        out = sys.stdout if args.output == "-" else open(args.output, "w")
        try:
            for root in roots:
                out.write(to_newick(root, decimals=args.decimals,
                                    label_internal=args.label_internal) + "\n")
        finally:
            if out is not sys.stdout:
                out.close()
        return

    # tskit format
    for i, root in enumerate(roots):
        ts = to_tree_sequence(root, sequence_length=args.sequence_length)
        if args.output == "-":
            sys.stdout.write(
                f"TreeSequence {i}: {ts.num_samples} samples, "
                f"{ts.num_nodes} nodes, {ts.num_edges} edges, "
                f"{ts.num_trees} tree(s), sequence_length={ts.sequence_length}\n")
        else:
            if args.replicates == 1:
                path = args.output
            else:
                base, dot, ext = args.output.rpartition(".")
                path = (f"{base}.{i}.{ext}" if dot else f"{args.output}.{i}")
            ts.dump(path)
            sys.stdout.write(f"wrote {path}\n")


if __name__ == "__main__":
    main()
