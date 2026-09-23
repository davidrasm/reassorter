#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run likelihood ratio test (LRT) for statistical support of reassortment events on simulated data
"""

from Espalier import MAF
from Espalier.RAxML import RAxMLRunner
from Espalier.SPRTester import SPRTester
import dendropy
import pandas as pd
import os
import glob
import sys
import datetime

# Utils.py is a standalone module (not inside the Espalier package) in the
# reference script — adjust this path to wherever your copy of it actually lives.
sys.path.insert(0, "/work/btb44/Espalier_testing/Espalier")  # <-- confirm this is where Utils.py is
import Utils

"Plot during debugging"
plot = True

"File names for input trees and alignments"
path = '/work/btb44/Espalier_testing/simulated_seqs_mu_0.01'
seg0_files = sorted(glob.glob(os.path.join(path, "arg_rate_*_seg0.fasta")))
print(f"Found {len(seg0_files)} seg0 files")

for seg0_seq_file in seg0_files:
    base = os.path.basename(seg0_seq_file).replace("_seg0.fasta", "")
    seg1_seq_file = seg0_seq_file.replace("_seg0.fasta", "_seg1.fasta")

    if not os.path.exists(seg1_seq_file):
        print(f"  Skipping {base}: no matching seg1 file")
        continue

    # True simulated trees (ground truth) — used only to pull tip labels/dates
    # for lsd dating, not to inform the ML search or its topology.
    true_seg0_tree_file = os.path.join(path, f"{base}_seg0.tre")
    true_seg1_tree_file = os.path.join(path, f"{base}_seg1.tre")

    if not os.path.exists(true_seg0_tree_file):
        print(f"  Skipping {base}: no true tree file found at {true_seg0_tree_file}")
        continue

    "Set up temp dir for temp tree output"
    temp_dir = 'temp-{date:%Y-%m-%d_%H%M%S}-{base}/'.format(
        date=datetime.datetime.now(), base=base)
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)

    "Initialize callable instances of Espalier objects"
    raxml = RAxMLRunner(raxml_path='raxml-ng', lsd_path='lsd', temp_dir=temp_dir)
    tester = SPRTester(raxml, temp_dir=temp_dir)

    "Tip dates and clock rate for lsd dating"
    tip_date_file = os.path.join(temp_dir, 'dates-lsd.txt')
    rate_file = os.path.join(temp_dir, 'rate.txt')
    mut_rate = 0.01  # match whatever rate was used to simulate this ARG set
    Utils.write_tip_dates(true_seg0_tree_file, tip_date_file)
    Utils.write_rate_file(mut_rate, rate_file)

    "Infer dated ML trees"
    seg0_tree_file = os.path.join(temp_dir, f"{base}_seg0.tre")
    seg1_tree_file = os.path.join(temp_dir, f"{base}_seg1.tre")
    raxml.get_dated_raxml_tree(seg0_seq_file, seg0_tree_file, tip_date_file, rate_file)
    raxml.get_dated_raxml_tree(seg1_seq_file, seg1_tree_file, tip_date_file, rate_file)

    "Load in trees"
    taxa = dendropy.TaxonNamespace()
    seg0_tree = dendropy.Tree.get(file=open(seg0_tree_file, 'r'), schema="newick",
                                   rooting="default-rooted", taxon_namespace=taxa)
    seg1_tree = dendropy.Tree.get(file=open(seg1_tree_file, 'r'), schema="newick",
                                   rooting="default-rooted", taxon_namespace=taxa)

    "Run LRT for seg0 and seg1"
    maf = MAF.get_maf_4cut(seg0_tree, seg1_tree, plot=False)
    if plot:
        print("MAF:")
        MAF.plot_maf(maf)

    results_df = tester(seg0_tree, seg0_seq_file, seg1_tree, seg1_seq_file, maf)
    results_df.to_csv(f"{base}_lrt_results.csv", index=False)
    print(f"  wrote {base}_lrt_results.csv")