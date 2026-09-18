#!/bin/bash

# ============================================================
# Parameters for MLE
# ============================================================

# Path to the MLE script
SCRIPT="/Users/blairblakeney/GitHub/reassorter/reassortment_MLE/20260914_BTB_reassortment_MLE.py"

# Directory where ARG simulation outputs are stored
ARG_DIR="/Users/blairblakeney/GitHub/reassorter/reassortment_MLE/arg_dataset"

# Name of output file containing true rates and MLEs
RESULTS_FILE="mle_summary.csv"

# Parameters passed to the MLE script (should match the simulations)
N_SAMPLES=10        
SEGMENTS=3          
NE=100.0            

# Optimizer bounds
LOWER_BOUND=0       
UPPER_BOUND=1       


# ============================================================
# Run MLE on each simulated ARG
# ============================================================

# Write header of output file
echo "true_rate,mle_rate" > "${RESULTS_FILE}"

for nodes_file in "${ARG_DIR}"/arg_rate_*.trees.nodes.csv
do
    # Base name shared by the nodes/edges files
    file_name=$(basename "${nodes_file}" .trees.nodes.csv)

    # True rate parsed from the filename
    true_rate=$(echo "${file_name}" | sed -E 's/arg_rate_//; s/\.trees$//')

    edges_file="${ARG_DIR}/${file_name}.trees.edges.csv"

    mle=$(python3 "${SCRIPT}" \
        -n "${N_SAMPLES}" \
        --segments "${SEGMENTS}" \
        -Ne "${NE}" \
        --lower-bound "${LOWER_BOUND}" \
        --upper-bound "${UPPER_BOUND}" \
        --nodes-csv "${nodes_file}" \
        --edges-csv "${edges_file}" \
        --reassortment-rate "${true_rate}")

    echo "${true_rate},${mle}" >> "${RESULTS_FILE}"
done
