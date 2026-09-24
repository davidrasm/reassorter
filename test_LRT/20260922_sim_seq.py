import glob
import os
import sys
import tskit

# Set path to wherever the Espalier GitHub was cloned
sys.path.insert(0, "/work/btb44/Espalier_testing/Espalier")
from Espalier.sim.ARGSimulator import sim_seqs

# Directory where simulated ARGs are stored
tree_dir = "/work/btb44/Espalier_testing/args_2_segments_mu_0.01"
tree_pattern = "*.trees"
output_dir = "/work/btb44/Espalier_testing/simulated_seqs_mu_0.01"

# Set parameters for simulating sequence data for each segment
segment_seq_length = 1000   # bp length per segment
mut_rate = 0.01
freqs = [0.25, 0.25, 0.25, 0.25]
kappa = 2.75

os.makedirs(output_dir, exist_ok=True)

trees_files = sorted(glob.glob(os.path.join(tree_dir, tree_pattern)))
print(f"Found {len(trees_files)} .trees files")

# Testing with just one tree
#trees_path = trees_files[0]

for trees_path in trees_files:
    base = os.path.splitext(os.path.basename(trees_path))[0]

    ts = tskit.load(trees_path)
    n_segments = ts.num_trees

    # Loop through all segment trees 
    for tr_num, tree in enumerate(ts.trees()):
        # A tree may span multiple segments if no reassortment event separated them.
        # But, seg0 and seg1 always get their own sequence files even when their trees are identical
        seg_start, seg_end = int(tree.interval[0]), int(tree.interval[1])

        for seg_idx in range(seg_start, seg_end):
            tree_file = os.path.join(output_dir, f"{base}_seg{seg_idx}.tre")
            seq_file = os.path.join(output_dir, f"{base}_seg{seg_idx}.fasta")

            # True tree for this segment
            with open(tree_file, "w") as f:
                print(tree.newick(), file=f)

            # Simulate an alignment along the true segment tree
            try:
                sim_seqs(
                    tree_file,
                    seq_file,
                    mut_rate=mut_rate,
                    seq_length=segment_seq_length,
                    freqs=freqs,
                    kappa=kappa,
                )
            except (Exception, SystemExit) as e:
                print(f"  FAILED on {tree_file}: {e}")