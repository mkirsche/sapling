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

# Runtimes
qlengths = [11, 21, 31, 41, 51, 101]
sapling01Runtimes = [49.4984, 135.788, 158.525, 162.19, 162.996, 166.227]
sapling1Runtimes = [22.0612, 106.697, 113.088, 116.37, 114.941, 114.533]
sapling25Runtimes = [25.8727, 81.8653, 91.242, 94.8738, 96.5893, 100.62]
runtimeNames = ['Bowtie', 'Mummer', 'Binary Search', 'Sapling (.01% overhead)', 'Sapling (1% overhead)', 'Sapling (25% overhead)']

bowtieTimes = [405, 558, 628, 737, 835, 1347]
mummerTimes = [870.457, 360.728, 220.997, 126.001, 109.334, 99.256]
bsTimes = [97.93, 288.3, 298.761, 306.561, 309.613, 299.606]

current_palette = sns.color_palette()

fig = plt.figure()
ax = plt.subplot(111)

sns.lineplot(qlengths, bowtieTimes, marker = "o", color = pal[0])
sns.lineplot(qlengths, mummerTimes, marker = "o", color = pal[1])
sns.lineplot(qlengths, bsTimes, marker = "o", color = pal[2])
sns.lineplot(qlengths, sapling01Runtimes, marker = "o", color = pal[3])
sns.lineplot(qlengths, sapling1Runtimes, marker = "o", color = pal[4])
sns.lineplot(qlengths, sapling25Runtimes, marker = "o", color = pal[5])
plt.ylim(-20, 1390)

plt.xlabel("Query Length")
plt.ylabel("$Runtime (seconds)$")
plt.title('Runtime with Variable Query Length')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.55, box.height])
plt.legend(labels = runtimeNames, bbox_to_anchor=(1.01, 0.92))
#plt.show()
plt.savefig('varqueryruntime.png')
