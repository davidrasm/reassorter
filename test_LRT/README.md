# Simulating ARGs with Reassortment and Estimating Reassortment Events with LRT

This workflow simulates ancestral recombination graphs (ARGs) with reassortment, simulates sequence data corresponding to these ARGs, generates maximum likelihood trees from this sequence data, and finally estimates reassortment events from these trees using a likelihood ratio test.

All scripts live in ~/GitHub/reassorter/test_LRT


## Overview

| Step | Task | Script |
|------|------|--------|
| 1 | Simulate ARGs | `20260903_BTB_arg.py`, `simulate_many_args.sh` |
| 2 | Clone Espalier dev-lrt branch|
| 3 | Simulate sequence data | `20260922_sim_seq.py` |
| 4 | Run likelihood ratio test on maximum likelihood trees | `run_sim_reassort_lrt.py` |
| 5 | Compare counts of true reassortment events versus estimated reassortment events | `compare_reassortment_counts.py` |

---

## 1. Simulate ARGs

	Note: In this case, I ran 100 ARG sims - with 2 segments, reassortment rates 0.001 to 1, Ne = 100, and 10 samples.


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

## 2. Clone Espalier dev-lrt branch

	- Run this from the command line while in your chosen conda environment

```bash
> git clone --branch dev-lrt https://github.com/davidrasm/Espalier.git
> conda activate py_env
> cd Espalier
> pip install -e .
```

---

## 3. Simulate sequence data

	Note: In this case, I used a mutation rate of 0.01.

	Script:** `20260922_sim_seq.py`

	- Make sure you have msprime and pyvolve installed in your conda environment

```bash
sbatch -t 1:00:00 --mem 4G -p common --wrap="python3 20260922_sim_seq.py"
```

---

## 4. Run likelihood ratio test on maximum likelihood trees

	Note: When running the below code, I got many warnings about 'Parent height is below sister subtree height' or 'Parent height is below subtree height.'

	Script: `run_sim_reassort_lrt.py`

	- Make sure you have RAxML loaded in your conda environment

```bash
sbatch -t 1:00:00 --mem 4G -p common --wrap="python3 run_sim_reassort_lrt.py"
```

---

## 5. Compare counts of true reassortment events versus estimated reassortment events

	Note: I used an alpha of 0.05 in this case.

	Script: `compare_reassortment_counts.py`

	Output: `reassortment_count_comparison.csv`

```bash
python3 compare_reassortment_counts.py
```

---

## Recommendations for next steps

A) Increase number of sites or increase mutation rate in order to reduce tree-inference error
B) Increase the number of samples from 10 so you don't hit a ceiling on the number of SPRs