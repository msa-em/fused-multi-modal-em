---
title: Multi-Modal Tutorial Dataset 2
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

:::{figure} ./figs/Figure_6_Output_2.png
:name: Raw vs Fused
:width: 700px
Comparison of raw input vs fused multi-modal Au-Fe$_3$O$_4$ HAADF elastic and EDX inelastic images
:::

```{warning} Step 0: Experimental Requirements
To reconstruct using fused multi-modal electron microscopy you need to collect both elastic (e.g. HAADF) and inelastic (e.g. EELS / EDX) maps of your material. For the elastic signal, it is important that it provides Z-contrast of your elements. For the inelastic signal, you should have all chemistries in your sample mapped. Solving for under-determined chemical maps, or using difficult to interpret elastic signals is outside the scope of this tutorial.
```

```{admonition} Step 1
Python Imports
```
:::{code-cell} ipython3
import fusion_utils as utils

from scipy.sparse import spdiags
import matplotlib.pyplot as plt
from tqdm import tqdm 
import numpy as np
import h5py
:::

```{admonition} Step 2
For this example, the dataset is stored in a .h5 file so this is how you can extract data and then save it numpy arrays for the fused multi-modal workflow.
```

:::{code-cell} ipython3
data = 'Au_Fe3O4_dataset.h5'

# Define element names and their atomic weights
elem_names=['Au', 'O', 'Fe']
elem_weights=[79,8,26]
# Parse elastic HAADF data and inelastic chemical maps based on element index from line above
with h5py.File(data, 'r') as h5_file:
    HAADF = np.array(h5_file['HAADF'])
xx = np.array([],dtype=np.float32)
for ee in elem_names:
	# Read chemical maps
    with h5py.File(data, 'r') as h5_file:
        chemMap = np.array(h5_file[ee])
        
  # Check if chemMap has the same dimensions as HAADF
	if chemMap.shape != HAADF.shape:
		raise ValueError(f"The dimensions of {ee} chemical map do not match HAADF dimensions.")
	

	# Set Noise Floor to Zero and Normalize Chemical Maps
	chemMap -= np.min(chemMap); chemMap /= np.max(chemMap)

	# Concatenate Chemical Map to Variable of Interest
	xx = np.concatenate([xx,chemMap.flatten()])
:::

```{admonition} Step 3
Reshape your data
```

```{danger} Caution!
Do not change the code below.
```

:::{code-cell} ipython3
# Make Copy of Raw Measurements for Poisson Maximum Likelihood Term 
xx0 = xx.copy()

gamma = 1.6 # Don't change this parameter from 1.6

# Image Dimensions
(nx, ny) = chemMap.shape; nPix = nx * ny
nz = len(elem_names)

# C++ TV Min Regularizers
reg = utils.tvlib(nx,ny)

# Data Subtraction and Normalization 
HAADF -= np.min(HAADF); HAADF /= np.max(HAADF)
HAADF=HAADF.flatten()

# Create Summation Matrix
A = utils.create_weighted_measurement_matrix(nx,ny,nz,elem_weights,gamma,1)
:::

```{admonition} Optional
:class: tip 
Plot your raw elastic/inelastic data
```

:::{code-cell} ipython3
fig, ax = plt.subplots(2, len(elem_names), figsize=(12, 8))
ax = ax.flatten()

for ii in range(len(elem_names)):
    ax[ii].imshow(xx0[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx, ny), cmap='gray')
    ax[ii].set_title(elem_names[ii])
    ax[ii].axis('off')
    
    ax[ii + len(elem_names)].imshow(xx0[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx, ny)[40:100, 50:110], cmap='gray')
    ax[ii + len(elem_names)].set_title(elem_names[ii] + ' Cropped')
    ax[ii + len(elem_names)].axis('off')

plt.show()
:::

```{admonition} Step 4
Fine tune your weights for each of the three parts of the cost function.
```

:::{code-cell} ipython3
# Convergence Parameters
lambdaHAADF = 1/nz # Do not modify this
lambdaChem = 0.4 # Ranges from 0 to 1, usually best around 0.2-0.3
nIter = 100 # Typically 10-15 will suffice
lambdaTV = 0.09; #Typically between 0.001 and 1
bkg = 0.55

# FGP TV Parameters
regularize = True; nIter_TV = 3; 
:::

```{admonition} Step 5
Run the Fused Multi-Modal algorithm
```

```{danger} Caution!
Do not change the code below.
```

:::{code-cell} ipython3
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
			xx[zz*nPix:(zz+1)*nPix] = reg.fgp_tv( xx[zz*nPix:(zz+1)*nPix].reshape(nx,ny), lambdaTV, nIter_TV).flatten()
			costTV[kk] += reg.tv( xx[zz*nPix:(zz+1)*nPix].reshape(nx,ny) )
			
	# Measure Cost Function
	costHAADF[kk] = lsqFun(xx); costChem[kk] = poissonFun(xx)
:::

```{admonition} Step 6
Assess Convergence
```

```{admonition} Be careful with $\lambda_{TV}$!
:class: attention
```

:::{code-cell} ipython3
# Display Cost Functions and Descent Parameters
utils.plot_convergence(costHAADF, lambdaHAADF, costChem, lambdaChem, costTV, lambdaTV)
# Show Reconstructed Signal
fig, ax = plt.subplots(2, len(elem_names), figsize=(12, 8))
ax = ax.flatten()

for ii in range(len(elem_names)):
    ax[ii].imshow(xx[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx, ny), cmap='gray')
    ax[ii].set_title(elem_names[ii])
    ax[ii].axis('off')
    
    ax[ii + len(elem_names)].imshow(xx[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx, ny)[40:100, 50:110], cmap='gray')
    ax[ii + len(elem_names)].set_title(elem_names[ii] + ' Cropped')
    ax[ii + len(elem_names)].axis('off')

plt.show()
:::

```{admonition} Step 7
Save your data
```

:::{code-cell} ipython3
save_folder_name='test'
utils.save_data(save_folder_name, xx0, xx, HAADF, A.dot(xx**gamma), elem_names, nx, ny, costHAADF, costChem, costTV, lambdaHAADF, lambdaChem, lambdaTV, gamma)
:::