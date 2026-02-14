import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

# ----------------- CLI arguments -----------------
p = argparse.ArgumentParser(
    description="Compare three runs (e.g. different rho values) on the same plot"
)

p.add_argument(
    "--run_dirs",
    nargs=3,
    required=True,
    help="Three directories, each containing sinkjet_history.csv",
)

p.add_argument(
    "--labels",
    nargs=3,
    required=True,
    help="Labels for the three runs (e.g. 1 2 3 or rho1 rho2 rho3)",
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

# ----------------- Helper function to read CSV -----------------
def load_history(run_dir):
    csv_path = os.path.join(run_dir, "sinkjet_history.csv")
    print("Reading:", csv_path)

    if not os.path.exists(csv_path):
        raise SystemExit(f"Could not find {csv_path}")

    cols = ["time", "Msink", "Mgas", "SFE_M0", "SFE_inst"]
    df = pd.read_csv(csv_path, comment="#", names=cols)

    # Athena can output duplicate times
    df = df.drop_duplicates(subset=["time"]).sort_values("time")

    return df

# ----------------- Load all runs -----------------
df1 = load_history(args.run_dirs[0])
df2 = load_history(args.run_dirs[1])
df3 = load_history(args.run_dirs[2])

# ----------------- Choose what to plot -----------------

if args.quantity == "Mtotal":
    y1 = df1["Msink"] + df1["Mgas"]
    y2 = df2["Msink"] + df2["Mgas"]
    y3 = df3["Msink"] + df3["Mgas"]
    ylabel = "Total Mass (Msink + Mgas)"
elif args.quantity in ["Msink", "Mgas", "SFE_M0", "SFE_inst"]:
    y1 = df1[args.quantity]
    y2 = df2[args.quantity]
    y3 = df3[args.quantity]
    ylabel = args.quantity
else:
    raise RuntimeError("Unknown quantity")

# ----------------- Plot -----------------
plt.figure()

plt.plot(
    df1["time"],
    y1,
    label=args.labels[0],
    linewidth=2,
)

plt.plot(
    df2["time"],
    y2,
    label=args.labels[1],
    linewidth=2,
    linestyle="--",
)

plt.plot(
    df3["time"],
    y3,
    label=args.labels[2],
    linewidth=2,
    linestyle=":",
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
    label1 = args.labels[0]
    label2 = args.labels[1]
    label3 = args.labels[2]
    args.outfile = f"{args.quantity}_compare_{label1}_{label2}_{label3}.png"

plt.savefig(args.outfile, dpi=150)
print(f"Wrote {args.outfile}")
