import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

# ----------------- CLI arguments -----------------
p = argparse.ArgumentParser(
    description="Compare multiple runs on the same plot"
)

p.add_argument(
    "--run_dirs",
    nargs="+",
    required=True,
    help="Directories containing sinkjet_history.csv",
)

p.add_argument(
    "--labels",
    nargs="+",
    required=True,
    help="Labels for each run (must match number of run_dirs)",
)

p.add_argument(
    "--quantity",
    choices=["SFE_M0", "SFE_inst", "Msink", "Mgas", "Mtotal"],
    required=True,
    help="Which quantity to plot",
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
    help="y-axis limits",
)

p.add_argument(
    "--outfile",
    default=None,
    help="Output PNG filename (default is auto-generated)",
)

args = p.parse_args()

# ----------------- Validation -----------------
if len(args.run_dirs) != len(args.labels):
    raise SystemExit("Number of run_dirs must match number of labels.")

# ----------------- Helper function -----------------
def load_history(run_dir):
    csv_path = os.path.join(run_dir, "sinkjet_history.csv")
    print("Reading:", csv_path)

    if not os.path.exists(csv_path):
        raise SystemExit(f"Could not find {csv_path}")

    cols = ["time", "Msink", "Mgas", "SFE_M0", "SFE_inst"]
    df = pd.read_csv(csv_path, comment="#", names=cols)

    df = df.drop_duplicates(subset=["time"]).sort_values("time")
    return df

# ----------------- Load all runs -----------------
dfs = [load_history(d) for d in args.run_dirs]

# ----------------- Choose quantity -----------------
if args.quantity == "Mtotal":
    ys = [df["Msink"] + df["Mgas"] for df in dfs]
    ylabel = "Total Mass (Msink + Mgas)"
elif args.quantity in ["Msink", "Mgas", "SFE_M0", "SFE_inst"]:
    ys = [df[args.quantity] for df in dfs]
    ylabel = args.quantity
else:
    raise RuntimeError("Unknown quantity")

# ----------------- Plot -----------------
plt.figure()

# Distinct styles
linestyles = ["-", "--", ":", "-."]
colors = plt.cm.tab20.colors  # 20 distinct colors

for i, (df, y, label) in enumerate(zip(dfs, ys, args.labels)):
    plt.plot(
        df["time"],
        y,
        label=label,
        linewidth=2,
        linestyle=linestyles[i % len(linestyles)],
        color=colors[i % len(colors)],
    )

plt.xlabel("Time [code units]")
plt.ylabel(ylabel)
plt.legend()
plt.tight_layout()

ax = plt.gca()
if args.xlim:
    ax.set_xlim(*args.xlim)
if args.ylim:
    ax.set_ylim(*args.ylim)

# ----------------- Save -----------------
if args.outfile is None:
    label_string = "_".join(args.labels)
    args.outfile = f"{args.quantity}_compare_{label_string}.png"

plt.savefig(args.outfile, dpi=150)
print(f"Wrote {args.outfile}")
