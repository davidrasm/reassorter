# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`reassorter` is a Python population genetics simulator. There are two standalone scripts, both outputting Newick or tskit TreeSequence format:

- **`coalescent.py`** — Kingman coalescent: simulates a single backward-time genealogy (one tree).
- **`arg.py`** — coalescent *with recombination*: simulates an ancestral recombination graph (ARG) and emits its local/marginal trees. Two interchangeable models via `--mode`: `hudson` (continuous genome, arbitrary breakpoints) and `reassortment` (segmented genome, whole-segment swapping, as in influenza).
- **`viz.py`** — matplotlib drawing of an ARG as a phylogenetic network (sample / coalescence / recombination nodes, time on the y-axis). Consumes the optional `ARGGraph` from `arg.py`; only imported when plotting.

`arg.py` started as a copy of `coalescent.py` and shares its conventions (Ne, ploidy, RNG, output formats), so keep terminology and CLI flags consistent across the two when editing.

## Running the Simulators

```bash
# Single coalescent tree
python3 coalescent.py -n 10 -Ne 1000 --ploidy 2 --format newick --seed 42 -r 5 -o out.tre

# ARG, Hudson recombination (default mode; -L required, --rho optional)
python3 arg.py -n 5 -Ne 1000 --rho 2e-6 -L 1000 --seed 42                 # local trees
python3 arg.py -n 5 -Ne 1000 --rho 2e-6 -L 1000 --format tskit -o arg.trees

# ARG, reassortment mode (--segments required)
python3 arg.py -n 5 -Ne 1000 --mode reassortment --segments 8 --reassortment-rate 1e-3
python3 arg.py -n 5 -Ne 1000 --mode reassortment --segments 8 --reassortment-rate 1e-3 --reassortment-bias 0.8

# Draw the ARG network to an image instead of writing trees (either mode)
# The --plot extension selects the format (png/pdf/svg/...; defaults to .png)
python3 arg.py -n 5 -Ne 1000 --rho 2e-6 -L 1000 --plot arg.png
python3 arg.py -n 5 -Ne 1000 --rho 2e-6 -L 1000 --plot arg.pdf

# Take defaults from a TOML config; CLI arguments still override it
python3 arg.py --config config.example.toml --seed 42
python3 arg.py --config config.example.toml -n 10 --rho 5e-6   # CLI wins
```

Both scripts take effective population size via `-Ne/--Ne`. Optional dependencies: `pip install tskit` (for `--format tskit`), `pip install matplotlib` (for `--plot`). `arg.py --config` parses TOML via the stdlib `tomllib`, so it requires Python 3.11+.

## Architecture

There is no build system, test suite, or packaging — each script is self-contained.

### coalescent.py

**Core data model**: `Node` (defined with `__slots__`) is the only data structure — `name`, `time` (generations before present), `children`. Leaf nodes are named `n<i>` (e.g. `n1`), internal nodes `i<j>`.

**Pipeline**:
```
simulate_coalescent(n, Ne, ploidy, rng)
    → root Node

root Node → to_newick(root)          # Newick string
         → to_tree_sequence(root)    # tskit TreeSequence (optional dep)
               └─ _tree_to_records(root)  # tskit-independent flatten step
```

- `simulate_coalescent`: backward-time loop — draws an exponential waiting time using rate `C(k,2) / (ploidy * Ne)`, randomly merges two of the `k` active lineages into a new parent node, until one lineage (the MRCA root) remains.
- `_tree_to_records`: tskit-independent core — flattens the Node tree into raw node/edge records so table contents can be validated without tskit installed.

### arg.py

Unlike `coalescent.py`, the ARG is **not** a `Node` tree — it is built directly as flat tskit-style tables, because recombination means different genomic intervals have different trees.

**Core data model**:
- A *lineage* is a list of `(left, right, node_id)` ancestral segments — the genome intervals it is ancestral over, with the node id ancestral on each. Segments are tracked because recombination splits a lineage's intervals between two parents.
- `node_time` / `node_is_sample`: per-node-id lists (samples are ids `0..n-1`).
- `edges`: flat `(left, right, parent, child)` records — the succinct/tskit representation. **Recombination is not a node**; it shows up as a child connecting to two parents over disjoint intervals.
- `OverlapCounter`: piecewise-constant count of how many lineages cover each position; used to detect when an interval has reached its MRCA (coverage drops to 1) so it can be dropped instead of climbing further.

**Pipeline**:
```
simulate_arg(n, Ne, rho, L, ploidy, rng)
    → (node_time, node_is_sample, edges, breakpoints)

tables → marginal_trees(...)   # list of (left, right, newick) local trees
       → to_tree_sequence(...) # tskit TreeSequence (optional dep)
edges  → squash_edges(edges)   # merge same parent/child across adjacent intervals
```

- `simulate_arg`: backward-time loop with coalescence (rate `C(k,2)/(ploidy*Ne)`, via `_merge_lineages`, recording edges only where lineages overlap) competing against a mode-dependent split event, until every position has coalesced. The `mode` branch is the *only* model-specific code; lineage selection is shared via `_weighted_choice`.
  - `mode="hudson"`: per-lineage rate `rho * span`; `_split_lineage` cuts at one uniform breakpoint over a continuous genome `[0, L)`.
  - `mode="reassortment"`: per-lineage rate is a constant `reassortment_rate` (0 unless the lineage carries ≥2 segments); `_reassort_split` sends each integer-aligned unit segment to one of two parents independently with prob `reassortment_bias`. Genome `L` is the segment count `K`. Events where all segments land on one parent are model-faithful no-ops (skipped, not recorded).
- `marginal_trees`: derives local trees per genomic interval from the edge table without tskit; merges adjacent intervals with identical topology.

**Explicit ARG graph (`record_graph=True`)**: the default succinct/tskit `edges` output has *no* recombination nodes and discards event times, so it can't draw a network. Passing `record_graph=True` to `simulate_arg` additionally builds an `ARGGraph` (namedtuples `ARGNode`/`ARGEdge`) with explicit sample, coalescence, and recombination nodes. The loop tracks `lineage_node[i]` (the graph node at the bottom of `pool[i]`'s current upward stretch) in lockstep with `pool`; each `ARGEdge` carries the ancestral segments it transmits, so a future per-lineage segment overlay needs no new bookkeeping. Node degrees: coalescence = 2 children/1 parent, recombination = 1 child/2 parents. Graph node ids are a separate id space from the tskit node ids, leaving the succinct output untouched.

**`simplify_arg(graph)`**: returns a copy with degree-2 nodes (one distinct parent + one distinct child) suppressed, splicing their two edges into one, iterated to a fixpoint. This collapses redundant "bubbles" — a recombination immediately undone by its two recombinant lineages re-coalescing, which appears as a recombination and coalescence node joined by a parallel edge pair and contributes nothing to any genealogy. Used by `viz.draw_arg(..., simplify=True)` (default) and toggled by the `--plot ... --no-simplify` CLI flag.

### viz.py

`draw_arg(graph, ...)` renders the `ARGGraph` with matplotlib in a rectangular phylogram style: lineages are vertical segments with a horizontal jog to the ancestor at its time. `_layout` assigns y = node time and x from a **spanning tree** of the ARG (each node claimed by the first parent reaching it in a downward DFS; leaves get in-order slots, internal nodes the mean of their tree children). Ordering leaves by topology gives every subtree a contiguous x-interval, so backbone horizontals never cross a vertical; only reticulation edges (a recombination node's second parent) can cross, which is unavoidable for an ARG. The two parental lineages leaving a recombination node fork left/right (by `fork`, the left parent to the left) so they read as two ancestors. Circles for sample/coalescence, squares for recombination annotated with breakpoint (hudson) or segment routing (reassortment). matplotlib is imported lazily. Legible only for small ARGs.

### Key Conventions (both scripts)

- Time increases going backward (samples at `t=0`)
- `Ne` is in individuals; diploid coalescence rate uses `ploidy * Ne` as the effective haploid size
- `rho` (arg.py hudson mode) is the per-site, per-generation recombination rate; a lineage's recombination rate scales with the span of ancestral material it carries
- `reassortment_bias` p and 1−p are equivalent (parents are exchangeable); p only controls how lopsided splits are — p=0.5 maximizes realized reassortment, p→0/1 makes events mostly no-ops
- Branch lengths = difference between parent and child `time`
- RNG is passed explicitly (`rng` parameter) for reproducibility
- `main` in each script is the argparse CLI entry point: replicates loop, output file management, optional tskit table building

### arg.py `--config` (TOML defaults)

`arg.py main` does a **two-phase parse**: `parse_known_args` pulls out `--config`, `_apply_config` loads the TOML and folds it in via `parser.set_defaults`, then `parse_args` runs. Precedence is `add_argument default < config file < CLI`. Because `set_defaults` cannot satisfy `required=True`, `-n/--num-samples` and `-Ne/--Ne` are *not* required at the argparse level (default `None`) and are checked for presence after parsing (so config can supply them). `set_defaults` also bypasses argparse's `type`/`choices` validation, so `_apply_config` validates keys against the parser's known dests and values against each action's `choices`. Config keys are dest names (e.g. `num_samples`, `genome_length`). See `config.example.toml`.
