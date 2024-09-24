import data.fusion_utils as utils
from scipy.sparse import spdiags
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm 
import numpy as np
import ipywidgets as widgets
from IPython.display import display

def return_reconstruction_plots(
    xx0,
    HAADF,
    A,
    bkg,
    n,
    elem_names,
    crop_roi,
    lambdaChem,
    lambdaTV,
    nIter,
    nIter_TV,
    regularize=True,
    gamma = 1.6,
    subtract_bkg = None,
):
    """ """

    xx = xx0.copy()
    if subtract_bkg is not None:
        xx = np.where((xx < subtract_bkg), 0, xx)
    
    nx, ny, nz = n
    sx0, sx1, sy0, sy1 = crop_roi
    nPix = nx * ny
    reg = utils.tvlib(nx,ny)
    lambdaHAADF = 1/nz
    
    costHAADF = np.zeros(nIter, dtype=np.float32)
    costChem = np.zeros(nIter, dtype=np.float32)
    costTV = np.zeros(nIter, dtype=np.float32)
    
    for kk in tqdm(range(nIter)):
        # Solve for the first two optimization functions $\Psi_1$ and $\Psi_2$
        xx -= gamma * spdiags(xx ** (gamma - 1), [0], nz * nx * ny, nz * nx * ny) * lambdaHAADF * A.transpose() * (A.dot(xx ** gamma) - HAADF) + lambdaChem * (1 - xx0 / (xx + bkg))

        # Enforce positivity constraint
        xx[xx < 0] = 0

        # FGP Regularization if turned on
        if regularize:
            for zz in range(nz):
                xx[zz * nPix:(zz + 1) * nPix] = reg.fgp_tv(xx[zz * nPix:(zz + 1) * nPix].reshape(nx, ny), lambdaTV, nIter_TV).flatten()
                # Measure TV Cost Function
                costTV[kk] += reg.tv(xx[zz * nPix:(zz + 1) * nPix].reshape(nx, ny))

        # Auxiliary Functions for measuring the cost functions
        lsqFun = lambda inData: 0.5 * np.linalg.norm(A.dot(inData ** gamma) - HAADF) ** 2
        poissonFun = lambda inData: np.sum(xx0 * np.log(inData + 1e-8) - inData)
        
        # Measure $\Psi_1$ and $\Psi_2$ Cost Functions
        costHAADF[kk] = lsqFun(xx)
        costChem[kk] = poissonFun(xx)

    utils.plot_convergence(costHAADF, lambdaHAADF, costChem, lambdaChem, costTV, lambdaTV)

    # Show Reconstructed Signal
    fig, ax = plt.subplots(2, len(elem_names) + 1, figsize=(12, 6.5))
    ax = ax.flatten()
    ax[0].imshow((A.dot(xx ** gamma)).reshape(nx, ny), cmap='gray')
    ax[0].set_title('HAADF')
    ax[0].axis('off')
    ax[1 + len(elem_names)].imshow((A.dot(xx ** gamma)).reshape(nx, ny)[70:130, 25:85], cmap='gray')
    ax[1 + len(elem_names)].set_title('HAADF Cropped')
    ax[1 + len(elem_names)].axis('off')

    for ii in range(len(elem_names)):
        ax[ii + 1].imshow(xx[ii * (nx * ny):(ii + 1) * (nx * ny)].reshape(nx, ny), cmap='gray')
        ax[ii + 1].set_title(elem_names[ii])
        ax[ii + 1].axis('off')
        ax[ii + 2 + len(elem_names)].imshow(xx[ii * (nx * ny):(ii + 1) * (nx * ny)].reshape(nx, ny)[sx0:sx1, sy0:sy1], cmap='gray')
        ax[ii + 2 + len(elem_names)].set_title(elem_names[ii] + ' Cropped')
        ax[ii + 2 + len(elem_names)].axis('off')

    fig.tight_layout()
    
    return None