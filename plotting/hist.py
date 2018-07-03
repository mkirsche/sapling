import sys
import matplotlib.pyplot as plt
import numpy as np

fn = sys.argv[1]

with open(fn) as f:
    vals = []
    lines = f.readlines()
    for l in lines:
        vals.append(int(l))
    borders = np.arange(0, 400, 2)
    plt.hist(vals, bins = borders)
    plt.title('Maximum overprediction error per segment')
    plt.show()
