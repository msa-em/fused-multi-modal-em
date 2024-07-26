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

```{admonition} Step 1
Python Imports
```
:::{code-cell} ipython3
from scipy.sparse import spdiags
import matplotlib.pyplot as plt
import fusion_utils as utils
from tqdm import tqdm 
import numpy as np
:::

```{admonition} Step 2: Load your data
```

:::{code-cell} ipython3
data = np.load('ZnSCu_dataset.npz')
# Define element names and their atomic weights
elem_names=['Cu', 'S', 'Zn', 'O']
elem_weights=[29,16,30, 8]
# Parse elastic HAADF data and inelastic chemical maps based on element index from line above
HAADF = data['HAADF']
xx = np.array([],dtype=np.float32)
for ee in elem_names:

	# Read Chemical Map for Element "ee"
	edsMap = data[ee]

	# Set Noise Floor to Zero and Normalize Chemical Maps
	edsMap -= np.min(edsMap); edsMap /= np.max(edsMap)

	# Concatenate Chemical Map to Variable of Interest
	xx = np.concatenate([xx,edsMap.flatten()])
:::

```{admonition} Step 3: Reshape your data
```

```{danger} Danger!
Do not change the code below.
```

:::{code-cell} ipython3
# Make Copy of Raw Measurements for Poisson Maximum Likelihood Term 
xx0 = xx.copy()

gamma = 1.6 # Don't change this parameter from 1.6

# Image Dimensions
(nx, ny) = edsMap.shape; nPix = nx * ny
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
fig, ax = plt.subplots(2,len(elem_names)+1,figsize=(12,8))
ax = ax.flatten()
ax[0].imshow(HAADF.reshape(nx,ny),cmap='gray'); ax[0].set_title('HAADF'); ax[0].axis('off')
ax[1+len(elem_names)].imshow(HAADF.reshape(nx,ny)[40:100,50:110],cmap='gray'); ax[1+len(elem_names)].set_title('HAADF Cropped'); ax[1+len(elem_names)].axis('off')

for ii in range(len(elem_names)):
    ax[ii+1].imshow(xx0[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny),cmap='gray'); ax[ii+1].set_title(elem_names[ii]); ax[ii+1].axis('off')
    ax[ii+2+len(elem_names)].imshow(xx0[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny)[40:100,50:110],cmap='gray'); ax[ii+2+len(elem_names)].set_title(elem_names[ii]+' Cropped'); ax[ii+2+len(elem_names)].axis('off')
plt.show()
:::

```{admonition} Step 4: Fine tune your weights for each of the three parts of the cost function.
```

:::{code-cell} ipython3
# Convergence Parameters
lambdaHAADF = 1/nz # Do not modify this
lambdaChem = 0.02
lambdaTV = 0.004; #Typically between 0.001 and 1
nIter = 100 # Typically 10-15 will suffice
bkg = 0.4

# FGP TV Parameters
regularize = True; nIter_TV = 3; 
:::

```{admonition} Step 5: Run the Fused Multi-Modal algorithm
```

```{danger} Danger!
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

```{admonition} Step 6: Assess Convergence
```

```{admonition} Be careful with $\lambda_{TV}$!
:class: attention
```

:::{code-cell} ipython3
# Display Cost Functions and Descent Parameters
utils.plot_convergence(costHAADF, lambdaHAADF, costChem, lambdaChem, costTV, lambdaTV)
# Show Reconstructed Signal
fig, ax = plt.subplots(2,len(elem_names)+1,figsize=(12,8))
ax = ax.flatten()
ax[0].imshow((A.dot(xx**gamma)).reshape(nx,ny),cmap='gray'); ax[0].set_title('HAADF'); ax[0].axis('off')
ax[1+len(elem_names)].imshow((A.dot(xx**gamma)).reshape(nx,ny)[40:100,50:110],cmap='gray'); ax[1+len(elem_names)].set_title('HAADF Cropped'); ax[1+len(elem_names)].axis('off')

for ii in range(len(elem_names)):
    ax[ii+1].imshow(xx[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny),cmap='gray'); ax[ii+1].set_title(elem_names[ii]); ax[ii+1].axis('off')
    ax[ii+2+len(elem_names)].imshow(xx[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny)[40:100,50:110],cmap='gray'); ax[ii+2+len(elem_names)].set_title(elem_names[ii]+' Cropped'); ax[ii+2+len(elem_names)].axis('off')
plt.show()
:::

```{admonition} Step 7: Save your data
```

:::{code-cell} ipython3
save_folder_name='test'
utils.save_data(save_folder_name, xx0, xx, HAADF, A.dot(xx**gamma), elem_names, nx, ny, costHAADF, costChem, costTV, lambdaHAADF, lambdaChem, lambdaTV, gamma)
:::