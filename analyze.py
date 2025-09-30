import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV (written by problem_sinkjet.cpp)
df = pd.read_csv('sinkjet_history.csv', comment='#',
                 names=['time','Msink','Mgas','SFE_M0','SFE_inst'])

print(df.head())  # quick sanity check
# Plot both definitions
plt.figure()
plt.plot(df['time'], df['SFE_M0'],   label='SFE = Msink / M0 (initial total)')
plt.plot(df['time'], df['SFE_inst'], label='SFE_inst = Msink / (Msink + gas)')
plt.xlabel('time')
plt.ylabel('star formation efficiency')
plt.legend()
plt.tight_layout()
#ax = plt.gca()

# pick limits that cover BOTH runs
#ax.set_xlim(0.00, 0.05)          # same time range
#ax.set_ylim(0.038475, 0.038650)     # example; choose what fits your data
## optional: consistent ticks/format
#from matplotlib.ticker import MaxNLocator, FormatStrFormatter
#ax.yaxis.set_major_locator(MaxNLocator(6))
#ax.yaxis.set_major_formatter(FormatStrFormatter('%.5f'))
#

plt.savefig('sfe_vs_time.png', dpi=150)
print('Wrote sfe_vs_time.png')
