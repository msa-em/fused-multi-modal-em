---
title: Code Tutorial
numbering:
  enumerator: 4.%s
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
### Guided Computation of Fused Multi-Modal Electron Microscopy

```{admonition} Step 1: Python Imports
First we import common python packages alongside some EM specific packages for i/o capabilites. [fusion_utils](https://github.com/jtschwar/Multi-Modal-2D-Data-Fusion/blob/170fea3292da7e6390bfff7236610eb0c8077ff7/EDX/fusion_utils.py) is a class of functions specifically built for 2D Fused MM.
```
:::{code-cell} ipython3
from scipy.sparse import spdiags
import matplotlib.pyplot as plt
import data/fusion_utils as utils
from tqdm import tqdm 
import numpy as np
import h5py
import hyperspy.api as hs
from skimage.io import imsave
:::

```{admonition} Step 2: Locate relevant EM data
Load and display indices of Elemental Maps and HAADF from raw EM files
```

:::{code-cell} ipython3
fname = '/data/example_dataset_1.emd'
# Load Raw Data and Reshape
file = hs.load(fname)
for idx, f in enumerate(file):
    print(f"{idx}: {f}")
:::

```{attention} Checking and extracting imported data
Make sure to keep track of which indices above correspond to your elements of interest
```

:::{code-cell} ipython3
# Parse Chemical Maps based on element index from line above
elementList = [7,4,3]
# List Element Names and Atomic Weights
elem_names=['Sc', 'Dy', 'O']
elem_weights=[21,66,8]
# Pull HAADF
HAADF=np.array(file[8]).flatten()
# Pull Elemental Maps
xx = np.array([],dtype=np.float32)
for ee in elementList:

	# Read Chemical Map for Element "ee"
	elemMap = np.array(file[ee])

	# Set Noise Floor to Zero and Normalize Chemical Maps
	elemMap -= np.min(elemMap); elemMap /= np.max(elemMap)

	# Concatenate Chemical Map to Variable of Interest
	xx = np.concatenate([xx,elemMap.flatten()])
:::

```{admonition} Step 3: Reshape relevant EM data
Pull Elemental Maps and HAADF from raw EM files.  Visualize the raw EM projections and chemical maps.
```

```{danger} Danger!
Do not change the code below.
```

:::{code-cell} ipython3
# Make Copy of Raw Measurements for Poisson Maximum Likelihood Term 
xx0 = xx.copy()

gamma = 1.6 # Don't change this parameter from 1.6

# Image Dimensions
(nx, ny) = elemMap.shape; nPix = nx * ny
nz = len(elementList)

# C++ TV Min Regularizers
reg = utils.tvlib(nx,ny)

# Data Subtraction and Normalization 
HAADF -= np.min(HAADF)
HAADF /= np.max(HAADF)

# Create Summation Matrix
A = utils.create_weighted_measurement_matrix(nx,ny,nz,elem_weights,gamma,1)

# Show Raw Data
utils.plot_elemental_images(xx, HAADF, elem_names, nx, ny, 1,(len(elem_names)+1))
:::

```{admonition} Step 4: Fine-tune and run Fused-MM Algorithm
The cost function parameters must be fine tuned in order to ensure convergence to a reasonable solution. Below is an example with its parameters finely tuned and clear convergence for all 3 functions is seen within 12 iterations of the algorithm. 
LambdaHAADF should not change from 1/nz
LambdaChem ranges from 0 to 1, we have found it works best when set to around 0.2-0.3
nIter can be set to your discretion, but most datasets we have worked with converge in <20 iterations
bkg works best between 1e-2 to 1e-1
ng works best around 3
lambdaTV typically works best between 0.001 and 1
```

:::{code-cell} ipython3
# Convergence Parameters
lambdaHAADF = 1/nz # Do not modify this
lambdaChem = 0.08 # Ranges from 0 to 1, usually best around 0.05-0.3
nIter = 30 # Typically 10-15 will suffice
bkg = 2.4e-1

# TV Min Parameters
regularize = True; ng = 3; 
lambdaTV = 0.15; #Typically between 0.001 and 1

xx = xx0.copy()

# Auxiliary Functions
lsqFun = lambda inData : 0.5 * np.linalg.norm(A.dot(inData**gamma) - HAADF) **2
poissonFun = lambda inData : np.sum(xx0 * np.log(inData + 1e-8) - inData)

# Main Loop
costHAADF = np.zeros(nIter,dtype=np.float32); costChem = np.zeros(nIter, dtype=np.float32); costTV = np.zeros(nIter, dtype=np.float32);
for kk in tqdm(range(nIter)):
	# HAADF Update
	xx -=  gamma * spdiags(xx**(gamma - 1), [0], nz*nx*ny, nz*nx*ny) * lambdaHAADF * A.transpose() * (A.dot(xx**gamma) - HAADF) + lambdaChem * (1 - xx0 / (xx + bkg))
	xx[xx<0] = 0

	# Regularization 
	if regularize:
		for zz in range(nz):
			xx[zz*nPix:(zz+1)*nPix] = reg.fgp_tv( xx[zz*nPix:(zz+1)*nPix].reshape(nx,ny), lambdaTV, ng).flatten()
			
	# Measure Cost Function
	costHAADF[kk] = lsqFun(xx); costChem[kk] = poissonFun(xx)

	# Measure Isotropic TV 
	if regularize:
		for zz in range(nz):
			costTV[kk] += reg.tv( xx[zz*nPix:(zz+1)*nPix].reshape(nx,ny) )


# Show Reconstructed Signal
utils.plot_elemental_images(xx,A.dot(xx**gamma),elem_names,nx,ny,1,4)

# Display Cost Functions and Descent Parameters
utils.plot_convergence(costHAADF, lambdaHAADF, costChem, lambdaChem, costTV, lambdaTV)
:::

```{admonition} Step 5: Save your Fused-MM data
Here we save the data as .tif files and record the cost function data to a .h5 file
```

:::{code-cell} ipython3
imsave('results/output_haadf.tif', (A.dot(xx**gamma)).reshape(nx,ny))
for ii in range(len(elementList)):
    imsave('results/output_{}.tif'.format(elementList[ii]), xx[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny))

# Save Cost Functions Output and Parameters (Experiment)
file = h5py.File('results/output_convergence.h5', 'w')
file.create_dataset('costHAADF', data=costHAADF); file['lambdaHAADF'] = lambdaHAADF
file.create_dataset('costChem', data=costChem); file['lambdaChem'] = lambdaChem
file.create_dataset('costTV', data=costTV); file['lambdaTV'] = lambdaTV
file['gamma'] = gamma; file.close()
:::