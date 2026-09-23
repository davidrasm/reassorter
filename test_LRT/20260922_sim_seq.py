import glob
import os
import sys
import tskit

sys.path.insert(0, "/work/btb44/Espalier_testing/Espalier")
from Espalier.sim.ARGSimulator import sim_seqs

# ---- CONFIG ----
tree_dir = "/work/btb44/Espalier_testing/args_2_segments_mu_0.01"
tree_pattern = "*.trees"
output_dir = "/work/btb44/Espalier_testing/simulated_seqs_mu_0.01"

mode = "reassortment"       # "hudson" or "reassortment" — must match how the ARGs were simulated
segment_seq_length = 1000   # [reassortment only] real bp length per segment (unit interval = 1 in ts coords)

mut_rate = 0.01
freqs = [0.25, 0.25, 0.25, 0.25]
kappa = 2.75
min_seg_length = 1          # [hudson only] skip trees narrower than this many bp
# ----------------

os.makedirs(output_dir, exist_ok=True)

trees_files = sorted(glob.glob(os.path.join(tree_dir, tree_pattern)))
print(f"Found {len(trees_files)} .trees files")

# Testing with just one tree
#trees_path = trees_files[0]

for trees_path in trees_files:
    base = os.path.splitext(os.path.basename(trees_path))[0]

    ts = tskit.load(trees_path)
    n_segments = ts.num_trees
    print(f"{base}: {n_segments} local tree(s)")

    for tr_num, tree in enumerate(ts.trees()):
        if mode == "reassortment":
            # A tree may span multiple segments if no reassortment event separated
            # them (ts coordinates are abstract unit indices [i, i+1), not bp).
            seg_start, seg_end = int(tree.interval[0]), int(tree.interval[1])

            with open("/tmp/_tmp_tree.tre", "w") as f:
                print(tree.newick(), file=f)

            for seg_idx in range(seg_start, seg_end):
                tree_file = os.path.join(output_dir, f"{base}_seg{seg_idx}.tre")
                seq_file = os.path.join(output_dir, f"{base}_seg{seg_idx}.fasta")

                with open(tree_file, "w") as f:
                    print(tree.newick(), file=f)

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

        else:
            # hudson: tree.interval is already in real bp coordinates
            seg_length = round(tree.interval[1] - tree.interval[0])
            if seg_length < min_seg_length:
                print(f"  skipping segment {tr_num} (length {seg_length} < min)")
                continue

            tree_file = os.path.join(output_dir, f"{base}_tree{tr_num}.tre")
            seq_file = os.path.join(output_dir, f"{base}_tree{tr_num}.fasta")

            with open(tree_file, "w") as f:
                print(tree.newick(), file=f)

            try:
                sim_seqs(
                    tree_file,
                    seq_file,
                    mut_rate=mut_rate,
                    seq_length=seg_length,
                    freqs=freqs,
                    kappa=kappa,
                )
            except (Exception, SystemExit) as e:
                print(f"  FAILED on {tree_file}: {e}")

print("Done.")