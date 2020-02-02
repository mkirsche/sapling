"""
File to train models. Testing and plotting is done in test.py, and the preprocessing of SA is done in preprocess_SA.py.
This code uses PyTorch 1.1.0.
"""

import torch
import torch.utils.data as utils_data
from torch.autograd import Variable
import numpy as np
from torch import nn, optim
import matplotlib.pyplot as plt
import argparse
from collections import OrderedDict
import os
torch.multiprocessing.set_start_method("spawn")

"""
USAGE:
python3 fit.py -s <layer size> -l <number of layers> -c <convergence window> -x <total number of chunks> -y < this model's chunk> -d < location/to/save/model/ > -f < location/of/preprocessed/suffix/array/ > 
e.g. python3 fit.py -s 8 -l 2 -c 20 -x 10 -y 1 -d /Users/user/NN/ -f /Users/user/SuffixArray/ 

Arguments:
-s: The size of each layer in the model
-l: The number of layers in the model
-c: The convergence window to be used
-t: The tolerance used (if improvement > tolerance in the convergence window is not seen, training ends)

-x: The total number of models being used on this data
-y: Which chunk of the data this model is trained on

-d: The output directory where the models and plots will be saved (must end with '/')
-f: The folder containing the preprocessed suffix array (must end with '/')

-p: If the loss should be plotted (or just saved)
-v: Verbose flag - if all output should be printed
"""

### ARGUMENTS ####
parser = argparse.ArgumentParser()

# Model Parameters
parser.add_argument("-s", required=False, default = 20,help="Layer Size")
parser.add_argument("-l", required=False, default = 1, help="Hidden Layers")
parser.add_argument("-c", required=False, default = 10, help="Convergence Window")
parser.add_argument("-e", required=False, default = 200, help="Max Epochs")
parser.add_argument("-t", required=False, default = 0.1, help="Convergence Threshold")
parser.add_argument("-x", required=False, default = 1, help ="Number of chunks to split the data into.")
parser.add_argument("-y", required=False, default = 1, help ="Which chunk this model is for (1,..., number of chunks).")

# I/O
parser.add_argument("-d", required=True, help="Output Directory")
parser.add_argument("-f", required=True, help ="Folder containing the preprocessed suffix array.")

# Output options - currently having some issues with plotting on electro, so leave at 0.
parser.add_argument("-p", required=False, default = 0, help ="Whether to generate plots - 1 if yes, 0 if no.")
parser.add_argument("-v", required=False, default = 0, help ="Print full output - 1 if yes, 0 if no.")
args = parser.parse_args()


### INITIALIZE ###

# Information on layer size, number of layers, convergence and number of epochs
layerSize = int(args.s)
hiddenLayers = int(args.l)
convergenceWindow = int(args.c)
num_epochs = int(args.e)
convergenceThreshold = float(args.t)

# numChunks = total number of chunks, currentChunkNum = which chunk this model is training on.
numChunks = int(args.x)
currentChunkNum = int(args.y)

# I/O
plotDir = str(args.d)
dataDir = str(args.f)

# Output options
plotFlag = int(args.p)
verbose = int(args.v)
if verbose == 0:
	v = False
else:
	v = True

# Print summary of the net
print("Using a net with " + str(hiddenLayers) + " hidden layers of size " + str(layerSize) + " with a window of " + str(convergenceWindow) + ", a threshold of " + str(convergenceThreshold) + " and " +str(num_epochs) + " max epochs.")

# Print the filename that will be used when saving this model
outputName = "s"+str(layerSize)+"_l"+str(hiddenLayers)+"_c"+str(convergenceWindow)+"_e"+str(num_epochs)+"_t"+str(convergenceThreshold)+"_x"+str(numChunks)+"_y"+str(currentChunkNum)
print("Model name: " + str(outputName))

# The directory for storing the models, plots and generated files
specificPlotDir = plotDir + outputName + "/"
if not os.path.exists(specificPlotDir):
	os.mkdir(specificPlotDir)

# The name of the .pkl file which will contain the model
modelDir = specificPlotDir+outputName+'.pkl'

# Print the information about the chunks
chunkString = "Training model " + str(currentChunkNum) + " of " + str(numChunks)
print(chunkString)

### CUDA CHECK ###
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print("Using "+ str(device))

### DATA LOADING - Load the saved numpy arrays from input folder ###

# X is the scaled decimal representations of the k-mers
Xfile = dataDir + 'x.npy'
X = np.load(Xfile)

# Y is the scaled suffix array positions
Yfile = dataDir + 'y.npy'
Y = np.load(Yfile)

# Res contains the normalized residual values
resfile = dataDir + 'res.npy'
res = np.load(resfile)

# True_res contains the original residual values (in terms of suffix array rows)
true_resfile = dataDir + 'true_res.npy'
true_res = np.load(true_resfile)

# Store some of the important values from the res array, so that we can convert error back into rows
res_min = np.min(true_res)
res_ptp = np.ptp(true_res)
min_pos = np.argmin(true_res)


### CHUNKING - by default, will do nothing to the arrays ###
print("Chunking data...")

# Divide the input data into chunks, and get appropriate chunk for this model
totalSize = len(X)
chunkSize = int(totalSize/numChunks)

# Getting the appropriate chunk.
if currentChunkNum != numChunks:
	start = (currentChunkNum-1)*chunkSize
	end = (currentChunkNum)*chunkSize
	X_chunk = X[start:end]
	Y_chunk = Y[start:end]
	res_chunk =res[start:end]
	true_res_chunk = true_res[start:end]
	print(start, end, len(X_chunk), len(Y_chunk))
else:
	start = (currentChunkNum-1)*chunkSize
	end = totalSize
	print(start, end)
	X_chunk = X[start:end]
	Y_chunk = Y[start:end]
	res_chunk =res[start:end]
	true_res_chunk = true_res[start:end]

### CREATING THE DATASETS ###

print("Loading data...")

# Store original, unscaled version of each
X_chunk_old = X_chunk
res_chunk_old = res_chunk

# Rescale so the chunk is between 0 and 1
# Edge case - if ptp = 0, set it to 1.
# Only happens if all values are the same, at which point they are all 0.
if np.ptp(X_chunk) == 0.0:
	X_chunk = (X_chunk - np.min(X_chunk))
else:
	X_chunk = (X_chunk - np.min(X_chunk))/np.ptp(X_chunk)

res_chunk = (res_chunk - np.min(res_chunk))/np.ptp(res_chunk)

# Get the arrays the model will be trained on
X_train = X_chunk
Y_train = res_chunk


# Load X and res values in to the dataloader, in batches of 64
batchsize = 64
training_samples = utils_data.TensorDataset(torch.from_numpy(X_train).float().to(device), torch.from_numpy(Y_train).float().to(device))
data_loader = utils_data.DataLoader(training_samples, batch_size=batchsize, shuffle=True, num_workers=0, pin_memory=False)

### MODEL ARCHITECTURE ###
# Model is created by adding iteratively adding layers and ReLu functions
# PyTorch cannot create a model directly from a list, so we create it from an ordered dictionary
architecture = []
layerCount = 0
architecture.append((str(layerCount), torch.nn.Linear(1, layerSize)))
layerCount += 1
architecture.append((str(layerCount), torch.nn.ReLU()))
layerCount += 1

for l in range(hiddenLayers-1):
	architecture.append((str(layerCount), torch.nn.Linear(layerSize, layerSize)))
	layerCount += 1
	architecture.append((str(layerCount), torch.nn.ReLU()))
	layerCount += 1

architecture.append((str(layerCount), torch.nn.Linear(layerSize,1)))

model = nn.Sequential(OrderedDict(architecture))

# Send the model to the device being used (CPU vs GPU)
model = model.to(device)

# Output the model specification
print(model)


### MODEL PARAMETERS ###
criterion = nn.MSELoss()
criterion = criterion.to(device)
optimizer = optim.Adam(model.parameters())#, lr=1e-3, momentum=0.8)
model.train()

### MODEL TRAINING ### 

print("Beginning training...")

# Store the loss values through the epochs
loss_list = []

# Tracking the best model, its loss, and the epoch in which we found it
best_loss = float('inf')
best_model = 0
best_epoch = 0

for epoch in range(num_epochs):

	# For each batch from the data, we will make a prediction, and compute the loss, and update the parameters
	for batch_idx, (data, target) in enumerate(data_loader):
		data, target = data.to(device), target.to(device)
		data, target = Variable(data), Variable(target)
		optimizer.zero_grad()
		output = model(data)
		output = output.to(device)
		loss = criterion(output, target.float())
		loss.backward()
		optimizer.step()

		# At the start of each epoch, print the current loss
		if batch_idx == 0:
			print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
			epoch, batch_idx * len(data), len(data_loader.dataset),
			100. * batch_idx / len(data_loader), loss.data.item()))

	# At the end of each epoch, print the final loss, and store the loss value
	print(loss.data.item())
	loss_list.append(loss.data.item())

	# Update the best model/loss/epoch
	if loss.data.item() < best_loss:
		best_loss = loss.data.item()
		best_model = model
		best_epoch = epoch
		torch.save(best_model, modelDir)

	# Checking for convergence - if we see no improvement in the convergenceWindow, we will break
	# Improvement is defined by the convergence threshold - by default, we look for 10% improvement every 50 epochs
	if epoch > convergenceWindow:
		# The best loss before this convergence window
		priorMin = np.min(loss_list[:(epoch-convergenceWindow)])
		# The best loss in this window
		windowMin = np.min(loss_list[(epoch-convergenceWindow):epoch])

		# The improvement observed in this window
		windowDiff = priorMin - windowMin

		# Print an update
		print ("Percent change: " +  str((windowDiff/priorMin)*100) + " from a priorMin of " + str(priorMin) + " to the current windowMin of " + str(windowMin))
		
		# If the improvment if not enough, we break.
		if windowDiff < convergenceThreshold*priorMin:
			print("Breaking at epoch " + str(epoch) + " with windowDiff " + str(windowDiff) + " and windowMin " + str(windowMin) + " and priorMin " + str(priorMin))
			break

# Save the loss numpy array to file.
loss_list = np.asarray(loss_list)
lossName = specificPlotDir + outputName + 'loss'
np.save(lossName, loss_list)

# Save the model
torch.save(best_model, modelDir)

# Switch into evaluating mode.
model.eval()

### SUMMARY - write Parameters and Details to file ###

# Human-readable output file
outFile_labels = specificPlotDir+outputName+'_labels.txt'
f = open(outFile_labels, "w")
f.write("Layer Size: " + str(layerSize) + '\n')
f.write("Hidden Layers: " + str(hiddenLayers)+ '\n')
f.write("Epochs: " + str(epoch)+ '\n')
f.write("Final Loss: "+ str(loss_list[-1])+ '\n')
f.write("Best Epoch: " + str(best_epoch)+ '\n')
f.write("Best Loss: " + str(best_loss)+ '\n')
f.close()

# File containing just the values
outFile = specificPlotDir+outputName+'.txt'
f = open(outFile, "w")
f.write(str(layerSize) + '\n')
f.write(str(hiddenLayers)+ '\n')
f.write(str(epoch)+ '\n')
f.write(str(loss_list[-1])+ '\n')
f.write(str(best_epoch)+ '\n')
f.write(str(best_loss)+ '\n')
f.close()

# Output
print("Layer Size: " + str(layerSize))
print("Hidden Layers: " + str(hiddenLayers))
print("Epochs: " + str(epoch))
print("Final Loss: "+ str(loss_list[-1]))
print("Best Epoch: " + str(best_epoch))
print("Best Loss: " + str(best_loss))

# Plotting epoch-by-epoch loss
if plotFlag == 1:
	# Plot the loss progression
	plt.plot(loss_list)
	plt.title("Loss")
	lossPlt = specificPlotDir+outputName+'_loss.png'
	# print(lossPlt)
	plt.savefig(lossPlt, bbox_inches='tight')
	# plt.show()
	plt.clf()