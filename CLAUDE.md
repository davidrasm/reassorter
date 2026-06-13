# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`reassorter` is a Python population genetics simulator. There are two standalone scripts, both outputting Newick or tskit TreeSequence format:

- **`coalescent.py`** — Kingman coalescent: simulates a single backward-time genealogy (one tree).
- **`arg.py`** — coalescent *with recombination* (Hudson's algorithm): simulates an ancestral recombination graph (ARG) and emits its local/marginal trees.

`arg.py` started as a copy of `coalescent.py` and shares its conventions (Ne, ploidy, RNG, output formats), so keep terminology and CLI flags consistent across the two when editing.

## Running the Simulators

```bash
# Single coalescent tree
python3 coalescent.py -n 10 -Ne 1000 --ploidy 2 --format newick --seed 42 -r 5 -o out.tre

# ARG with recombination (--rho and -L are required)
python3 arg.py -n 5 -Ne 1000 --rho 2e-6 -L 1000 --seed 42                 # local trees
python3 arg.py -n 5 -Ne 1000 --rho 2e-6 -L 1000 --format tskit -o arg.trees
```

Both scripts take effective population size via `-Ne/--Ne`. Optional tskit output requires: `pip install tskit`

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

- `simulate_arg`: backward-time loop with two competing events — coalescence (rate `C(k,2)/(ploidy*Ne)`, via `_merge_lineages`, recording edges only where lineages overlap) and recombination (per-lineage rate `rho * span`, via `_split_lineage` at a uniform breakpoint). Runs until every position has coalesced.
- `marginal_trees`: derives local trees per genomic interval from the edge table without tskit; merges adjacent intervals with identical topology.

### Key Conventions (both scripts)

- Time increases going backward (samples at `t=0`)
- `Ne` is in individuals; diploid coalescence rate uses `ploidy * Ne` as the effective haploid size
- `rho` (arg.py) is the per-site, per-generation recombination rate; a lineage's recombination rate scales with the span of ancestral material it carries
- Branch lengths = difference between parent and child `time`
- RNG is passed explicitly (`rng` parameter) for reproducibility
- `main` in each script is the argparse CLI entry point: replicates loop, output file management, optional tskit table building
