import pandas as pd
import matplotlib.pyplot as plt

#df = pd.read_csv("mle_summary.csv")
df = pd.read_csv("mle_summary_start_0.5.csv")

fig, ax = plt.subplots(figsize=(6, 6))

ax.scatter(df["true_rate"], df["mle_rate"], alpha=0.6, edgecolor="k", linewidth=0.3)

# 1:1 reference line (perfect recovery)
lims = [
    min(df["true_rate"].min(), df["mle_rate"].min()),
    max(df["true_rate"].max(), df["mle_rate"].max()),
]
ax.plot(lims, lims, "r--", linewidth=1, label="y = x")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("True reassortment rate")
ax.set_ylabel("MLE reassortment rate")
ax.set_title("MLE recovery of reassortment rate")
ax.legend()
ax.set_aspect("equal", adjustable="box")

plt.tight_layout()
plt.savefig("true_vs_mle_rate_start_0.5.png", dpi=300)
plt.show()