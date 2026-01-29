import argparse, os
import pandas as pd
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--run_dir", default=".", help="folder containing sinkjet_history.csv")
p.add_argument("--xlim", type=float, nargs=2, default=None, help="e.g. --xlim 0 0.5")
p.add_argument("--ylim", type=float, nargs=2, default=None, help="e.g. --ylim -1e-4 2e-4")
p.add_argument("--outfile", default="mass_change_vs_time.png")
p.add_argument("--percent", action="store_true",
               help="plot percent change instead of absolute change")
args = p.parse_args()

csv_path = os.path.join(args.run_dir, "sinkjet_history.csv")
if not os.path.exists(csv_path):
    raise SystemExit(f"Could not find {csv_path}. Did the run write it?")

cols = ["time","Msink","Mgas","SFE_M0","SFE_inst"]
df = pd.read_csv(csv_path, comment="#", names=cols)
df = df.drop_duplicates(subset=["time"]).sort_values("time")

print("First few rows:\n", df.head())
print("Last few rows:\n", df.tail())

# --- compute changes relative to initial values ---
Ms0 = df["Msink"].iloc[0]
Mg0 = df["Mgas"].iloc[0]

if args.percent:
    dMs = 100.0 * (df["Msink"] - Ms0) / (Ms0 if Ms0 != 0 else 1.0)
    dMg = 100.0 * (df["Mgas"]  - Mg0) / (Mg0 if Mg0 != 0 else 1.0)
    dMt = 100.0 * ((df["Msink"] + df["Mgas"]) - (Ms0 + Mg0)) / ((Ms0 + Mg0) if (Ms0 + Mg0) != 0 else 1.0)
    ylab = "Δ mass [%]"
else:
    dMs = df["Msink"] - Ms0
    dMg = df["Mgas"]  - Mg0
    dMt = (df["Msink"] + df["Mgas"]) - (Ms0 + Mg0)
    ylab = "Δ mass [code units]"

# --- plot ---
plt.figure()
plt.plot(df["time"], dMs, label="ΔMsink")
plt.plot(df["time"], dMg, label="ΔMgas", linestyle="--")
plt.plot(df["time"], dMt, label="ΔMtotal", linestyle=":")
plt.xlabel("Simulation time (code units)")
plt.ylabel(ylab)
plt.legend()
plt.tight_layout()

ax = plt.gca()
if args.xlim: ax.set_xlim(*args.xlim)
if args.ylim:
    ax.set_ylim(*args.ylim)
else:
    # give a small padding so near-flat lines are visible
    vmin = min(dMs.min(), dMg.min())
    vmax = max(dMs.max(), dMg.max())
    if vmin == vmax:
        pad = 1e-6 if not args.percent else 1e-3
        ax.set_ylim(vmin - pad, vmax + pad)

out_path = os.path.join(args.run_dir, args.outfile)
plt.savefig(out_path, dpi=150)
print(f"Wrote {out_path}")

