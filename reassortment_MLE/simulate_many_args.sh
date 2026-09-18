#!/bin/bash

# ============================================================
# Parameters for simulations
# ============================================================

# Path to the simulation script
SCRIPT="/Users/blairblakeney/GitHub/reassorter/20260903_BTB_arg.py"

# Reassortment rates: draw N_RATES values log-uniformly between LOW and HIGH
N_RATES=100
LOW=0.001
HIGH=1

# Simulation parameters passed to the script
N_SAMPLES=10        # -n
NE=100              # -Ne
MODE="reassortment" # --mode
SEGMENTS=3          # --segments
FORMAT="tskit"      # --format


# ============================================================
# Draw N_RATES reassortment rates from a log-uniform distribution
# ============================================================
reassortment_rates=()
while IFS= read -r rate; do
    reassortment_rates+=("$rate")
done < <(python3 -c "
import numpy as np
np.random.seed()
log_low, log_high = np.log10(${LOW}), np.log10(${HIGH})
rates = 10 ** np.random.uniform(log_low, log_high, size=${N_RATES})
for r in rates:
    print(f'{r:.6f}')
")

# ============================================================
# Loop through rates and run/submit a job for each
# ============================================================
for rate in "${reassortment_rates[@]}"
do
    # Slurm version 
    # sbatch -t 00:10:00 \
    #     --job-name="rate_${rate}" \
    #     --output="${rate}.out" \
    #     --error="${rate}.err" \
    #     --wrap="python3 ${SCRIPT} -n ${N_SAMPLES} -Ne ${NE} --mode ${MODE} --segments ${SEGMENTS} --reassortment-rate ${rate} --format ${FORMAT} -o arg_rate_${rate}.trees --plot arg_rate_${rate}.png"

    # Local version
    python3 "${SCRIPT}" \
        -n "${N_SAMPLES}" \
        -Ne "${NE}" \
        --mode "${MODE}" \
        --segments "${SEGMENTS}" \
        --reassortment-rate "${rate}" \
        --format "${FORMAT}" \
        -o "arg_rate_${rate}.trees" \
        --plot "arg_rate_${rate}.png"
done
