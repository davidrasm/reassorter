import re
import ast
import glob
import pandas as pd
import matplotlib.pyplot as plt

def classify_recombination(meta_str):
    meta = ast.literal_eval(meta_str)
    to_first = meta.get('to_first', [])
    to_second = meta.get('to_second', [])
    return 'observable' if (to_first and to_second) else 'unobservable'

records = []

# NOTE: adjust this glob + regex to match however your script actually
# names its nodes.csv files (e.g. "arg_rate_0.0123_nodes.csv")
for nodes_file in glob.glob("arg_rate_*.trees.nodes.csv"): 
    match = re.search(r"arg_rate_([\d.eE+-]+).trees.nodes\.csv", nodes_file)
    rate = float(match.group(1))

    df = pd.read_csv(nodes_file)  # use sep="\t" instead if these are tab-delimited
    recomb = df[df['type'] == 'recombination']

    counts = recomb['meta'].apply(classify_recombination).value_counts()
    n_obs = counts.get('observable', 0)
    n_unobs = counts.get('unobservable', 0)

    records.append({
        'rate': rate,
        'observable': n_obs,
        'unobservable': n_unobs,
        'total': n_obs + n_unobs,
    })

results = pd.DataFrame(records).sort_values('rate')
results.to_csv('reassortment_summary.csv', index=False)

fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(results['rate'], results['total'], label='Total', alpha=0.7)
ax.scatter(results['rate'], results['observable'], label='Observable', alpha=0.7)
ax.scatter(results['rate'], results['unobservable'], label='Unobservable', alpha=0.7)
ax.set_xscale('log')
ax.set_xlabel('Assigned reassortment rate (log scale)')
ax.set_ylabel('Number of reassortment events')
ax.legend()
plt.tight_layout()
plt.savefig('reassortment_rate_vs_events.png', dpi=300)