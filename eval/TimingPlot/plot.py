'''
A script for plotting runtime of Sapling with different parameters across many genome lengths
Requires hard-coding the values since each one comes from its own experiment
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math

sns.set()

# Genome lengths
xs = [4641652, 100286401, 230481012, 782475302, 2934876451, 14271578887]
logxs = [math.log(i, 2) for i in xs]

# Runtimes
sapling00001 = [64.1779, 114.843, 116.667, 131.352, 155.422, 222.989]
saplingOnePercent = [24.7783, 72.6055, 73.4093, 77.2295, 113.859, 179.87]
#saplingTenPercent = [19.5705,53.525, 60.0477, 60.2346, 91.2241, 123.117]
saplingTwentyfivePercent = [19.007, 48.6434, 57.5504, 58.5787, 83.7056, 115.208]
binSearch = [93.8439, 173.31, 192.317, 222.014, 288.3, 369.465]
runtimeNames = ['Binary Search', 'Sapling 0.01% Overhead', 'Sapling 1% Overhead', 'Sapling 25% Overhead']

sns.lineplot(logxs, binSearch, marker = "o")
sns.lineplot(logxs, sapling00001, marker = "o")
sns.lineplot(x = logxs, y = saplingOnePercent, marker = "o")
sns.lineplot(logxs, saplingTwentyfivePercent, marker = "o")

plt.xlabel("$log_{2}(genome size)$")
plt.ylabel("$Runtime (seconds)$")
plt.title('Runtime Scaling with Genome Size')
plt.legend(labels = runtimeNames)
plt.savefig('genomelengthruntime.png')
