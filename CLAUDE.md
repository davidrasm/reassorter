# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`reassorter` is a single-file Python population genetics simulator implementing the Kingman coalescent model. It simulates backward-time genealogies and outputs phylogenetic gene trees in Newick or tskit TreeSequence format.

## Running the Simulator

```bash
# Basic usage (required flags only)
python3 coalescent.py -n <num_samples> -N <effective_size>

# Common options
python3 coalescent.py -n 10 -N 1000 --ploidy 2 --format newick --seed 42 -r 5 -o out.tre
```

Optional tskit output requires: `pip install tskit`

## Architecture

Everything lives in `coalescent.py`. There is no build system, test suite, or packaging.

### Core Data Model

`Node` (defined with `__slots__`) is the only data structure. Each node carries `name`, `time` (generations before present), and `children`. Leaf nodes are named `n<i>` (e.g. `n1`), internal nodes `i<j>`.

### Simulation Pipeline

```
simulate_coalescent(n, Ne, ploidy, rng)
    → root Node

root Node → to_newick(root)          # Newick string
         → to_tree_sequence(root)    # tskit TreeSequence (optional dep)
               └─ _tree_to_records(root)  # tskit-independent flatten step
```

**`simulate_coalescent`**: Backward-time loop — at each step, draws an exponential waiting time using rate `C(k,2) / (ploidy * Ne)`, randomly merges two of the `k` active lineages into a new parent node at the new time, until one lineage (the MRCA root) remains.

**`_tree_to_records`**: Separates tskit integration concerns — flattens the Node tree into raw node/edge record lists without importing tskit, so it can be tested independently.

**`main`**: argparse-based CLI entry point. Handles replicates loop, output file management, and optional tskit table building.

### Key Conventions

- Time increases going backward (leaves at `t=0`, root at largest `t`)
- `Ne` is in individuals; diploid coalescence rate uses `ploidy * Ne` as the effective haploid size
- Branch lengths in Newick = difference between parent and child `time` values
- RNG is passed explicitly (`rng` parameter) to enable reproducible simulation
