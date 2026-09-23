#!/usr/bin/env python3
"""
Compare reassortment counts per simulated ARG:
  - observable / unobservable / total events counted directly from the ARG nodes.csv
  - SPRs in the MAF of the TRUE segment trees (what a topology-based method could detect)
  - SPRs in the MAF of the inferred ML trees, and how many the LRT supports
"""

import ast
import glob
import os
import pandas as pd
import dendropy
from Espalier import MAF

alpha = 0.05
results_dir = '/work/btb44/Espalier_testing/lrt_test_results_mu_0.01'  # folder containing the *_lrt_results.csv files
true_dir = '/work/btb44/Espalier_testing/simulated_seqs_mu_0.01'
nodes_dir = '/work/btb44/Espalier_testing/args_2_segments_mu_0.01'


def classify_recombination(meta_str):
    """Observable if the two segments go to different parents, unobservable otherwise."""
    meta = ast.literal_eval(meta_str)
    to_first = meta.get('to_first', [])
    to_second = meta.get('to_second', [])
    return 'observable' if (to_first and to_second) else 'unobservable'


def count_arg_events(base):
    """Count reassortment events in the true ARG for this simulation."""
    nodes_file = os.path.join(nodes_dir, f'{base}.trees.nodes.csv')
    if not os.path.isfile(nodes_file):
        return None, None
    nodes = pd.read_csv(nodes_file)
    recomb = nodes[nodes['type'] == 'recombination']
    counts = recomb['meta'].apply(classify_recombination).value_counts()
    return int(counts.get('observable', 0)), int(counts.get('unobservable', 0))

rows = []
for csv in sorted(glob.glob(os.path.join(results_dir, '*_lrt_results.csv'))):
    base = os.path.basename(csv).replace('_lrt_results.csv', '')
    df = pd.read_csv(csv)

    # Estimated: SPRs in the MAF of the ML trees, and how many the LRT supports
    n_spr_ml = len(df)
    n_sig = int((df['p_value'] < alpha).sum())
    n_sig_bonf = int((df['p_value'] < alpha / max(n_spr_ml, 1)).sum())

    # "Truth" comparable to the method: SPRs in the MAF of the true segment trees
    taxa = dendropy.TaxonNamespace()
    t0 = dendropy.Tree.get(path=os.path.join(true_dir, f'{base}_seg0.tre'), schema='newick',
                           rooting='default-rooted', taxon_namespace=taxa)
    t1 = dendropy.Tree.get(path=os.path.join(true_dir, f'{base}_seg1.tre'), schema='newick',
                           rooting='default-rooted', taxon_namespace=taxa)
    true_maf = MAF.get_maf_4cut(t0, t1, plot=False)
    n_spr_true = len(true_maf) - 1

    # Ground truth: events counted directly from the ARG
    n_obs, n_unobs = count_arg_events(base)
    if n_obs is None:
        print(f"  Warning: no nodes.csv found for {base}")

    rate = float(base.replace('arg_rate_', '')) if base.startswith('arg_rate_') else None
    rows.append({'base': base, 'reassort_rate': rate,
                 'ARG_observable': n_obs, 'ARG_unobservable': n_unobs,
                 'ARG_total': None if n_obs is None else n_obs + n_unobs,
                 'true_SPRs': n_spr_true, 'ML_SPRs': n_spr_ml,
                 'LRT_supported': n_sig, 'LRT_supported_bonf': n_sig_bonf})

out = pd.DataFrame(rows).sort_values('reassort_rate')
out['error_vs_trueSPRs'] = out['LRT_supported'] - out['true_SPRs']
out['error_vs_ARG_observable'] = out['LRT_supported'] - out['ARG_observable']
out.to_csv('reassortment_count_comparison.csv', index=False)

print(out.to_string(index=False))
print(f"\nAlignments compared: {len(out)}")
for truth, err in [('true_SPRs', 'error_vs_trueSPRs'),
                   ('ARG_observable', 'error_vs_ARG_observable')]:
    print(f"\nLRT_supported vs {truth}:")
    print(f"  Exact matches: {(out[err] == 0).mean():.2%}")
    print(f"  Mean error (estimated - true): {out[err].mean():.2f}")
    print(f"  Correlation: {out[truth].corr(out['LRT_supported']):.3f}")
