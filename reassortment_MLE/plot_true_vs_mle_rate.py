import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# ============================================================
# Command-line arguments
# ============================================================
parser = argparse.ArgumentParser(
    description="Compare true vs. MLE reassortment rates: print error/correlation "
                "statistics and save a scatter plot."
)
parser.add_argument("csv", help="Input CSV with 'true_rate' and 'mle_rate' columns")
parser.add_argument(
    "-o", "--output",
    help="Output plot filename (default: same name as the CSV, with a .png extension)",
)
args = parser.parse_args()

output = args.output or str(Path(args.csv).with_suffix(".png"))

# Load data with true rate versus MLE rate
df = pd.read_csv(args.csv)

x = df["true_rate"]
y = df["mle_rate"]

# Drop rows with missing values
mask = x.notna() & y.notna()
x, y = x[mask], y[mask]

# Statistics
# Mean squared error (raw scale)
mse = np.mean((y - x) ** 2)
rmse = np.sqrt(mse)

# Pearson correlation (raw scale)
r, p = stats.pearsonr(x, y)

# Spearman rank correlation 
rho, p_rho = stats.spearmanr(x, y)

print(f"MSE: {mse:.4g}  RMSE: {rmse:.4g}")
print(f"Pearson r: {r:.4f} (p={p:.3g})")
print(f"Spearman rho: {rho:.4f} (p={p_rho:.3g})")

# Log-scale statistics (only rows where both values are > 0)
pos = (x > 0) & (y > 0)
lx, ly = np.log10(x[pos]), np.log10(y[pos])

log_mse = np.mean((ly - lx) ** 2)
log_r, _ = stats.pearsonr(lx, ly)

print(f"Log10 MSE: {log_mse:.4g}")
print(f"Log10 Pearson r: {log_r:.4f}")
print(f"Rows dropped from log-scale stats: {(~pos).sum()}")

# Plot true vs. estimated rate
fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(x, y, alpha=0.6, edgecolor="k", linewidth=0.3)

# 1:1 reference line (perfect recovery)
lims = [min(x.min(), y.min()), max(x.max(), y.max())]
ax.plot(lims, lims, "r--", linewidth=1, label="y = x")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("True reassortment rate")
ax.set_ylabel("MLE reassortment rate")
ax.set_title("MLE recovery of reassortment rate")
ax.legend()
ax.set_aspect("equal", adjustable="box")

plt.tight_layout()
plt.savefig(output, dpi=300)
print(f"Saved plot to {output}")