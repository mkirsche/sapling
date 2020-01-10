"""
Plots a sample of the suffix array for each of six genome
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import random
import math
import sys
import seaborn as sns

ofn = sys.argv[1]

if ofn == 'gc.png':
  fns = ['suffixArraySim1.txt', 'suffixArraySim2.txt', 'suffixArraySim3.txt', 'suffixArraySim4.txt', 'suffixArraySim5.txt', 'suffixArraySim6.txt']
  names = ['Uniform Random', '60% GC-Content', '75% GC-Content', '90% GC-Content', '100% GC-Content', 'All C\'s']
else:
  fns = ['suffixArraySim7.txt', 'suffixArraySim8.txt', 'suffixArraySim9.txt', 'suffixArraySim10.txt', 'suffixArraySim11.txt', 'suffixArraySim12.txt']
  names = ['CT Repeat', 'CAT Repeat', 'ACGT Repeat', 'ACTTCA Repeat', 'Length 16 Repeat', 'Length 50 Repeat']
  ofn = 'repeats.png'

sns.set()
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(28,14))
plt.subplots_adjust(wspace=0.1, hspace=0.15)
plt.suptitle('Suffix Array Distributions', fontsize=16, y = .94, fontweight='heavy')
minys = []
maxys = []
for fileindex in range(0, len(fns)):
  fn = fns[fileindex]
  name = names[fileindex]
  xs = []
  ys = []
  with open(fn) as f:
    for line in f:
      tokens = line.split()
      xs.append(int(tokens[1]))
      ys.append(int(tokens[2]))
  minys.append(min(ys))
  maxys.append(max(ys))
  plot = sns.scatterplot(xs, ys, c = ['rebeccapurple' for i in range(0, len(xs))], ax = axs[fileindex/3][fileindex%3], linewidth=0)
  if fileindex / 3 == 1:
    plot.set_xlabel('k-mer', fontsize=14)
  if fileindex%3 == 0:
    plot.set_ylabel('Suffix Array Position', fontsize=14)
  plot.set_title(name, fontsize=16)
  #axs[fileindex/3][fileindex%3].set_yticklabels([])
  #axs[fileindex/3][fileindex%3].set_xticklabels([])
  axs[fileindex/3][fileindex%3].set_xlim(-.2e12, 5e12)
  #plt.clf()
  #plt.cla()
  #plt.close()

fig.savefig(ofn)

