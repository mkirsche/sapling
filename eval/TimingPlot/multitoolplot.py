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
saplingRuntimes = [291.782, 277.642 , 261.552, 251.606, 244.442, 252.177, 251.159, 190.957, 209.15, 183.053,
203.781, 167.8, 168.62, 159.528, 152.803, 131.041, 119.709, 142.572, 114.502, 133.928, 137.158, 114.152, 
135.96, 87.2305]
runtimeNames = ['Bowtie', 'Mummer', 'Binary Search', 'Sapling']

bowtieTime = 564
mummerTime = 396.905
binarySearchTime = 301.056

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
