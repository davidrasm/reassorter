# Simulating ARGs with Reassortment and Estimating Reassortment Rates with MLE

This workflow simulates ancestral recombination graphs (ARGs) with reassortment, then estimates the reassortment rate from each ARG by maximum likelihood.

All scripts live in ~/GitHub/reassorter/reassortment_MLE

The data and output files from simulating 100 ARGs and conducting MLE are stored in ~/GitHub/reassorter/reassortment_MLE/arg_dataset
  - These were the results for the relationship between MLEs and true reassortment rates:
    	MSE: 0.002031  RMSE: 0.04507
	    Pearson r: 0.9822 (p=5.73e-73)
	    Spearman rho: 0.9861 (p=3.79e-78) 


## Overview

| Step | Task | Script |
|------|------|--------|
| 1 | Simulate ARGs | `20260903_BTB_arg.py`, `simulate_many_args.sh` |
| 2 | Plot observed reassortment events vs. input rate | `reassort_counts_versus_rates.py` |
| 3 | Estimate rates by maximum likelihood | `20260914_BTB_reassortment_MLE.py`, `mle_for_many_args.sh` |
| 4 | Plot MLEs vs. true rates; compute MSE and correlation | `plot_true_vs_mle_rate.py` |

---

## 1. Simulate ARGs

	Script: `20260903_BTB_arg.py`

	Single ARG simulation at the command line:

```bash
python3 20260903_BTB_arg.py \
    -n 10 -Ne 100 \
    --mode reassortment --segments 3 \
    --reassortment-rate 1e-2 \
    --format tskit \
    -o arg_w_tables.trees \
    --plot arg_w_tables.png
```

	Many ARGs simulated at once using bash script:

	- Simulation parameters (number of samples, Ne, mode, number of segments, reassortment rate range, and output format) are set at the top of the bash script.

```bash
bash ~/GitHub/reassorter/reassortment_MLE/simulate_many_args.sh
```

---

## 2. Plot Observed Reassortment Events vs. Input Rate

	Script:** `reassort_counts_versus_rates.py`

	- Run this from the directory containing the ARG `.csv` outputs.

```bash
python3 ~/GitHub/reassorter/reassortment_MLE/reassort_counts_versus_rates.py
```

---

## 3. Estimate Reassortment Rates (MLE)

	Script: `20260914_BTB_reassortment_MLE.py`

	Single ARG estimation at the command line:

```bash
python3 ~/GitHub/reassorter/reassortment_MLE/20260914_BTB_reassortment_MLE.py \
    -n 10 --segments 3 -Ne 100.0 \
    --lower-bound 0 --upper-bound 1 \
    --nodes-csv arg_w_tables.trees.nodes.csv \
    --edges-csv arg_w_tables.trees.edges.csv \
    --reassortment-rate 0.01
```


	MLE for many ARGs using bash script:

	- MLE parameters are set at the top of the bash script. The output is a `.csv` with `true_rate` and `mle_rate` columns.
	- Currently set up so that the starting value for the optimizer is the true rate used in the ARG simulation. This doesn’t have to be the case for MLE to work!


```bash
bash ~/GitHub/reassorter/reassortment_MLE/mle_for_many_args.sh
```

---

## 4. Plot MLEs vs. True Rates

	Script: `plot_true_vs_mle_rate.py`

	- Takes the MLE `.csv` from step 3 as the first argument and the output plot filename via `-o`. Prints MSE, RMSE, Pearson r, and Spearman rho (raw and log10 scale), then saves the plot.

```bash
python3 ~/GitHub/reassorter/reassortment_MLE/plot_true_vs_mle_rate.py \
    mle_summary.csv \
    -o true_vs_mle_rate.png
```

	
