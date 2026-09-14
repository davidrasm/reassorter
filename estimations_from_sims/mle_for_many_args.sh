#!/bin/bash

# Retrieve all simulated arg output files
arg_output_files=

for file in arg_output_files
do
    # get beginning of file name
    file_name= 

    sbatch -t 01:00:00 \ 
    -job-name=mle_rate_${rate} \
    --output=mle_${rate}.out \
    --error=mle_${rate}.err \
    --wrap="python3 20260914_BTB_reassortment_MLE.py -n 10 --segments 3 --reassortment-rate ${rate} -Ne 100.0 \
      --lower-bound 0 --upper-bound 1 \
      --nodes-csv ${file_name}.nodes.csv --edges-csv ${file_name}.edges.csv"

      # Write out MLE to a file! Maybe "true" rate in one column and MLE in the other
done
