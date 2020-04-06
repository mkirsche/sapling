'''
A script for plotting runtime of Sapling compared
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math

sns.set()

pal = sns.color_palette()

# Genome lengths
logxs = [i for i in range(6, 30)]

# Runtimes
saplingRuntimes = [272.504, 271.494, 283.534, 195.717, 213.098, 199.751, 193.502, 190.024, 178.008,
  208.957, 193.096, 193.263, 155.422, 152.319, 132.068, 162.492, 120.197, 116.291, 113.859, 106.215,
99.927, 89.7132, 91.2241, 83.7056]
runtimeNames = ['Bowtie', 'Mummer', 'Binary Search', 'Sapling']

bowtieTime = 558
mummerTime = 360.728
binarySearchTime = 288.3

current_palette = sns.color_palette()

fig = plt.figure()
ax = plt.subplot(111)
ax.axhline(y = bowtieTime, c = pal[0])
ax.axhline(y = mummerTime, c = pal[1])
ax.axhline(y = binarySearchTime, c = pal[2])
sns.lineplot(logxs, saplingRuntimes, marker = "o", color = pal[3])
plt.ylim(-20, 620)

plt.xlabel("$log_{2}(Number of Bins)$")
plt.ylabel("$Runtime (seconds)$")
plt.title('Runtime of Exact Matching')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])
plt.legend(labels = runtimeNames, bbox_to_anchor=(1.01, 0.92))
#plt.show()
plt.savefig('multitoolruntime.png')
