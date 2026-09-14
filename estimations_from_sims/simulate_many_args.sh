#!/bin/bash

import numpy
# Draw 100 reassortment rates from a uniform distribution on a log scale
reassortment_rates = log_uniform_sample() # from 0.001 to 1

# Simulate args

for rate in reassortment_rates
do
    sbatch -t 01:00:00 \ 
    -job-name=rate_${rate} \
    --output=${rate}.out \
    --error=${rate}.err \
    --wrap="python3 20260903_BTB_arg.py -n 10 -Ne 100 --mode reassortment --segments 3 --reassortment-rate ${rate} --format tskit -o arg_rate_${rate}.trees --plot arg_rate_${rate}.png"
done
