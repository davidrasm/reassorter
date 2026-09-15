#!/bin/bash

# Directory containing simulated ARG outputs (nodes/edges CSVs)
#ARG_DIR="/Users/blairblakeney/GitHub/reassorter/estimations_from_sims/arg_dataset"

# Retrieve all simulated ARG node files
#for nodes_file in "${ARG_DIR}"/arg_rate_*.trees.nodes.csv
#do
    # strip path and "_nodes.csv" suffix to get the shared base name
#    file_name=$(basename "${nodes_file}" .trees.nodes.csv)

    # extract the true rate from the filename, e.g. arg_rate_0.012345 -> 0.012345
#    rate=$(echo "${file_name}" | sed -E 's/arg_rate_//')

#    edges_file="${ARG_DIR}/${file_name}.trees.edges.csv"

    #sbatch -t 01:00:00 \
    #    --job-name=mle_rate_${rate} \
    #    --output=mle_${rate}.out \
    #    --error=mle_${rate}.err \
    #    --wrap="python3 20260914_BTB_reassortment_MLE.py -n 10 --segments 3 -Ne 100.0 \
    #      --lower-bound 0 --upper-bound 1 \
    #      --nodes-csv ${nodes_file} --edges-csv ${edges_file} \
    #      --reassortment-rate ${rate} \
#done


ARG_DIR="/Users/blairblakeney/GitHub/reassorter/estimations_from_sims/arg_dataset"
RESULTS_FILE="mle_summary_start_0.5.csv"

# write header once
echo "true_rate,mle_rate" > "${RESULTS_FILE}"

for nodes_file in "${ARG_DIR}"/arg_rate_*.trees.nodes.csv
do
    file_name=$(basename "${nodes_file}" .trees.nodes.csv)
    #rate=$(echo "${file_name}" | sed -E 's/arg_rate_//; s/\.trees$//')
    true_rate=$(echo "${file_name}" | sed -E 's/arg_rate_//; s/\.trees$//')
    rate=0.5

    edges_file="${ARG_DIR}/${file_name}.trees.edges.csv"

    # capture whatever the python script prints to stdout
    mle=$(python3 /Users/blairblakeney/GitHub/reassorter/20260914_BTB_reassortment_MLE.py -n 10 --segments 3 -Ne 100.0 \
        --lower-bound 0 --upper-bound 1 \
        --nodes-csv "${nodes_file}" --edges-csv "${edges_file}" --reassortment-rate "${rate}")

    echo "${true_rate},${mle}" >> "${RESULTS_FILE}"
done
