'''
A script for plotting runtime of Bowtie with different sampling frequencies
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
samplingfreqs = [1, 2, 4, 8, 16, 32, 64]
runtimes = [315, 335, 336, 354, 412, 530, 761]

sapling01time = 155.422
sapling1time= 113.859
sapling25time = 83.7056

runtimeNames = ['Bowtie', 'Sapling 25%', 'Sapling 1%', 'Sapling .01%']
current_palette = sns.color_palette()

fig = plt.figure()
ax = plt.subplot(111)

sns.lineplot(samplingfreqs, runtimes, marker = "o", color = pal[0])
ax.axhline(y = sapling25time, c = pal[1])
ax.axhline(y = sapling1time, c = pal[2])
ax.axhline(y = sapling01time, c = pal[3])

plt.ylim(-20, 1390)

plt.xlabel("Bowtie Sampling Frequency")
plt.ylabel("$Runtime (seconds)$")
plt.title('Bowtie Runtimes')

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.95, box.height])
plt.legend(labels = runtimeNames, bbox_to_anchor=(1.01, 0.92))
#plt.show()
plt.savefig('bowtiefreqruntime.png')
