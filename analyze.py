import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

# ----------------- CLI arguments -----------------
p = argparse.ArgumentParser()
p.add_argument(
    "--run_dir",
    default=".",
    help="folder containing sinkjet_history.csv",
)
p.add_argument(
    "--xlim",
    type=float,
    nargs=2,
    default=None,
    help="x-axis limits, e.g. --xlim 0 5",
)
p.add_argument(
    "--ylim",
    type=float,
    nargs=2,
    default=None,
    help="y-axis limits, e.g. --ylim 0.03 0.05",
)
p.add_argument(
    "--outfile",
    default="sfe_vs_time.png",
    help="output PNG filename (saved inside run_dir)",
)

args = p.parse_args()

# ----------------- Read CSV -----------------
csv_path = os.path.join(args.run_dir, "sinkjet_history.csv")
print("DEBUG: reading", csv_path)

if not os.path.exists(csv_path):
    raise SystemExit(f"Could not find {csv_path}")

cols = ["time","Msink","Mgas","SFE_M0","SFE_inst"]
df = pd.read_csv(csv_path, comment="#", names=cols)

# remove duplicate times (Athena can print same time multiple times)
df = df.drop_duplicates(subset=["time"]).sort_values("time")

print("First few rows:\n", df.head())
print("Last few rows:\n", df.tail())

print("\nSFE_M0: min = ", df["SFE_M0"].min(), "max = ", df["SFE_M0"].max())
print("SFE_inst: min = ", df["SFE_inst"].min(), "max = ", df["SFE_inst"].max())

if df["SFE_M0"].isna().all():
    print("NOOOO:SFE_M0 column is N/A")
if df["SFE_inst"].isna().all():
    print("NOOOO: SFE_inst column is N/A")
    
# ----------------- Plot SFE only -----------------
plt.figure()
plt.plot(df["time"], df["SFE_M0"],
         label="SFE_M0 = Msink / M0",
         linestyle="--", linewidth=2)

plt.plot(df["time"], df["SFE_inst"],
         label="SFE_inst = Msink / (Msink + Mgas)",
         linestyle="-", linewidth=2)

plt.xlabel("time [code units]")
plt.ylabel("Star Formation Efficiency")
plt.legend()
plt.tight_layout()

ax = plt.gca()
if args.xlim:
    ax.set_xlim(*args.xlim)
if args.ylim:
    ax.set_ylim(*args.ylim)

# ----------------- Save -----------------
out_path = os.path.join(args.run_dir, args.outfile)
os.makedirs(args.run_dir, exist_ok=True)
plt.savefig(out_path, dpi=150)
print(f"Wrote {out_path}")
