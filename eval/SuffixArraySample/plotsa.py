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

fns = ['ecoli.sa.sample', 'celegans.sa.sample', 'chr1.sa.sample', 'tomato.sa.sample', 'human.sa.sample', 'wheat.sa.sample']
names = ['E. coli', 'C. elegans', 'H. sapiens (hg38 chr1)', 'S. lycopersicum', 'H. sapiens (hg38)', 'T. aestivum']
shortnames = ['ecoli', 'celegans', 'chr1', 'tomato', 'human', 'wheat']
sns.set()
fig, axs = plt.subplots(ncols=3, nrows=2, figsize=(28,14))
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.suptitle('Suffix Array Distributions', fontsize=16, y = .94, fontweight='heavy')
minys = []
maxys = []
for fileindex in range(0, len(fns)):
  fn = fns[fileindex]
  name = names[fileindex]
  shortname = shortnames[fileindex]
  xs = []
  ys = []
  with open(fn) as f:
    for line in f:
      tokens = line.split()
      xs.append(int(tokens[1]))
      ys.append(int(tokens[2]))
  minys.append(min(ys))
  maxys.append(max(ys))
  plot = sns.lineplot(xs, ys, color = 'rebeccapurple', ax = axs[fileindex/3][fileindex%3])
  if fileindex / 3 == 1:
    plot.set_xlabel('k-mer', fontsize=14)
  if fileindex%3 == 0:
    plot.set_ylabel('Suffix Array Position', fontsize=14)
  plot.set_title(name, fontsize=16)
  axs[fileindex/3][fileindex%3].set_yticklabels([])
  axs[fileindex/3][fileindex%3].set_xticklabels([])
  #axs[fileindex/3][fileindex%3].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
  #axs[fileindex/3][fileindex%3].ticklabel_format(style='sci', scilimits=(0,0))
  #plot.get_figure().savefig(shortname + '_sa.png')
  #plt.clf()
  #plt.cla()
  #plt.close()


numticks = 10
for i in range(0, 6):
  axs[i/3][i%3].set_yticks(np.linspace(0, maxys[i], numticks))
  ticks = axs[i/3][i%3].get_yticks()
  axs[i/3][i%3].set_ylim(-.2 * ticks[1], 1.2 * ticks[-1] - .2 * ticks[-2])

fig.savefig('all_sa.png')

