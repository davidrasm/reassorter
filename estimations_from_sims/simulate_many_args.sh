#!/bin/bash

# --- Draw N reassortment rates from a log-uniform distribution ---
N=100
LOW=0.001
HIGH=1

reassortment_rates=()
while IFS= read -r rate; do
    reassortment_rates+=("$rate")
done < <(python3 -c "
import numpy as np
np.random.seed()
log_low, log_high = np.log10(${LOW}), np.log10(${HIGH})
rates = 10 ** np.random.uniform(log_low, log_high, size=${N})
for r in rates:
    print(f'{r:.6f}')
")

#mapfile -t reassortment_rates < <(python3 -c "
#import numpy as np
#np.random.seed()  # remove or set a fixed int for reproducibility
#log_low, log_high = np.log10(${LOW}), np.log10(${HIGH})
#rates = 10 ** np.random.uniform(log_low, log_high, size=${N})
#for r in rates:
#    print(f'{r:.6f}')
#")

# --- Loop through rates and submit a job for each ---
for rate in "${reassortment_rates[@]}"
do
    #sbatch -t 01:00:00 \
    #    --job-name=rate_${rate} \
    #    --output=${rate}.out \
    #    --error=${rate}.err \
    #    --wrap="python3 20260903_BTB_arg.py -n 10 -Ne 100 --mode reassortment --segments 3 --reassortment-rate ${rate} --format tskit -o arg_rate_${rate}.trees --plot arg_rate_${rate}.png"
    python3 /Users/blairblakeney/GitHub/reassorter/20260903_BTB_arg.py -n 10 -Ne 100 --mode reassortment --segments 3 --reassortment-rate ${rate} --format tskit -o arg_rate_${rate}.trees --plot arg_rate_${rate}.png
done
