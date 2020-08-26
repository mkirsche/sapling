'''
A script for plotting memory usage of different tools
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math

sns.set()

sns.set()
current_palette = sns.color_palette()
mems = [12.29, 68.69, 79.861, 80.119, 88.246]
names = ['Bowtie', 'Mummer', 'Sapling .01%', 'Sapling 1%', 'Sapling 25%']

sns.barplot(names, mems)

plt.xlabel("Software")
plt.ylabel("$Memory (GB)$")
plt.title('Memory Usage')

plt.savefig('memory.png')
