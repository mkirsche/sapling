import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
fn = sys.argv[1]
pal = sns.color_palette()
with open(fn) as f:
  toPlot = []
  names = []
  goodness = []
  xs = []
  ys = []
  ps = []
  sns.set()
  for line in f:
    tokens = line.split(' ')
    if len(tokens) == 1:
      numBins = int(tokens[0])
      for i in range(0, numBins):
        toPlot.append([])
        xs.append([])
        ys.append([])
        ps.append([])
        names.append('')
        goodness.append(0)
    else:
      binId = int(tokens[0])
      plotNum = int(tokens[1])
      val = int(tokens[2])
      xs[plotNum].append(int(tokens[4]))
      ys[plotNum].append(int(tokens[5]))
      ps[plotNum].append(int(tokens[6]))
      toPlot[plotNum].append(val)
      names[plotNum] = str(binId)
      goodness[plotNum] = int(tokens[3])
  for i in range(0, len(toPlot)):
    clr = pal[2]
    #sns.distplot(toPlot[i], kde=False, bins = 50, color=clr)
    #plt.title('bin ' + names[i])
    #plt.savefig('figures/binHist' + str(i+1) + '.png')
    #plt.cla()
    #plt.clf()
    #plt.close()
    sns.lineplot(x=xs[i], y=ys[i], color = pal[0])
    sns.lineplot(x=xs[i], y=ps[i], color = clr)
    plt.title('bin ' + names[i])
    plt.savefig('figures/binScatter' + str(i+1) + '.png')
    plt.cla()
    plt.clf()
    plt.close()

