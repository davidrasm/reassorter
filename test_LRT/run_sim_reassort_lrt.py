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
import datetime
import glob

"Plot during dubugging"
plot = True

"File names for input trees and alignments"
path = '/work/btb44/Espalier_testing/simulated_seqs_mu_0.01'
seg0_files = glob.glob(os.path.join(path, "arg_rate_*_seg0.fasta"))
for seg0_seq_file in seg0_files:
    base = os.path.basename(seg0_seq_file).replace("_seg0.fasta", "")
    seg1_seq_file = seg0_seq_file.replace("_seg0.fasta", "_seg1.fasta") # but with seg1 instead of seg0


    "Set up temp dir for temp tree output"
    temp_dir = 'temp-{date:%Y-%m-%d_%H%M%S}-{base}/'.format(
        date=datetime.datetime.now(), base=base)
    if not os.path.isdir(temp_dir):
        os.mkdir(temp_dir)
        
    "Initialize callable instances of Espalier objects"
    raxml = RAxMLRunner(raxml_path='raxml-ng-mpi',lsd_path='lsd',temp_dir=temp_dir)
    tester = SPRTester(raxml,temp_dir=temp_dir)

    "Infer ML trees"
    seg0_tree_file = os.path.join(temp_dir, f"{base}_seg0.tre")
    seg1_tree_file = os.path.join(temp_dir, f"{base}_seg1.tre")
    raxml.get_raxml_tree(seg0_seq_file, seg0_tree_file)
    raxml.get_raxml_tree(seg1_seq_file, seg1_tree_file)

    "Load in trees"
    taxa = dendropy.TaxonNamespace()
    seg0_tree = dendropy.Tree.get(file=open(seg0_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)
    seg1_tree = dendropy.Tree.get(file=open(seg1_tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa)

    """
        Run LRT for seq0 and seg1
    """
    maf = MAF.get_maf_4cut(seg0_tree, seg1_tree, plot=False)
    if plot:
        print("MAF:")
        MAF.plot_maf(maf)
    results_df = tester(seg0_tree,seg0_seq_file,seg1_tree,seg1_seq_file,maf)
    results_df.to_csv(f"{base}_lrt_results.csv", index=False)
 