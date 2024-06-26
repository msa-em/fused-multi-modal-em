from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from skimage.io import imsave
from tqdm import tqdm
import numpy as np
import tifffile
import h5py

def save_data(folder_name, chem_original, chem_fused, haadf_original, haadf_fused, eList, nx, ny,costHAADF,costChem,costTV,lambdaHAADF,lambdaChem,lambdaTV,gamma):
    file = h5py.File(folder_name+'/Fused_Reconstruction.h5', 'w')
    file.create_dataset('raw_haadf', data=haadf_original)
    file.create_dataset('fused_haadf', data=haadf_fused)
    imsave(folder_name+'/haadf_recon.tif', haadf_fused.reshape(nx,ny))
    for ii in range(len(eList)):
        file.create_dataset('raw_'+eList[ii], data=(chem_original[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny)))
        file.create_dataset('fused_'+eList[ii], data=(chem_fused[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny)))
        imsave(folder_name+'/{}_signal.tif'.format(eList[ii]), chem_fused[ii*(nx*ny):(ii+1)*(nx*ny)].reshape(nx,ny))
    file.create_dataset('costHAADF', data=costHAADF); file['lambdaHAADF'] = lambdaHAADF
    file.create_dataset('costChem', data=costChem); file['lambdaChem'] = lambdaChem
    file.create_dataset('costTV', data=costTV); file['lambdaTV'] = lambdaTV
    file['gamma'] = gamma; file.close()

def plot_convergence(costLS, L_LS, costP, L_PS, costTV, L_TV):
	fig, (ax1, ax2, ax3) = plt.subplots(3,1,figsize=(12,6))

	nIter = costLS.shape[0]

	ax1.plot(costLS)
	ax1.set_ylabel(r'$\frac{1}{2} ||Ax^{\gamma} - b||^2$')
	ax1.set_title(r'$\mathbf{\lambda_{HAADF}: }$' + f'{L_LS:.3f}')  # Updated line
	ax1.set_xticklabels([])
	ax1.set_xlim([0, nIter-1])
	ax1.grid(True)


	ax2.plot(costP)
	ax2.set_ylabel(r'$\sum (x - b \cdot \log(x))$')
	ax2.set_title(r'$\mathbf{\lambda_{Chem}: }$' + f'{L_PS:.3f}')
	ax2.set_xticklabels([])
	ax2.set_xlim([0,nIter-1])
	ax2.grid(True)

	ax3.plot(costTV)
	ax3.set_ylabel(r'$\sum ||x||_{TV}$')
	ax3.set_xlabel('Iteration #')
	ax3.set_title(r'$\mathbf{\lambda_{TV}: }$' + f'{L_TV:.3f}')
	ax3.set_xlim([0,nIter-1])
	ax3.grid(True)

	plt.show()

def create_weighted_measurement_matrix(nx, ny, nz, zNums, gamma,method=0):
	#Create Measurement Matrix.
	vals = np.zeros([nz*ny*nx], dtype=np.float16)
	row =  np.zeros([nz*ny*nx], dtype=int)
	col =  np.zeros([nz*ny*nx], dtype=int)
	vals[:] = 1

	ii = 0; ind = 0
	while ii < nz*nx*ny:
		for jj in range(nz):
			row[ii+jj] = ind
			col[ii+jj] = ind + nx*ny*jj
	
			if method == 0:
				pass 
			if method == 1:
				vals[ii+jj] = zNums[jj] / np.mean(zNums) # Z_{i}^{gamma} / ( sum_{i} ( Z_i ) / Nelements )
			if method == 2:
				vals[ii+jj] = zNums[jj]**gamma / np.mean(np.power(zNums, gamma)) # Z_{i}^{gamma} / ( sum_{i} ( Z_i^{gamma} ) / Nelements )
			if method == 3:
				vals[ii+jj] = zNums[jj] / np.sum(zNums) # Z_{i} / sum_{i} Z_i
			if method == 4:
				vals[ii+jj] = zNums[jj]**gamma / np.sum(np.power(zNums, gamma)) # Z_{i}^{gamma} / sum_{i} ( Z_i^{gamma} )

		ii += nz
		ind += 1
	A = csr_matrix((vals, (row, col)), shape=(nx*ny, nz*nx*ny), dtype=np.float32)

	return A

# Class for Gradient Projection Method for Total Variation Denoising 
class tvlib:

    def __init__(self, nx, ny):
        self.nx = nx; self.ny = ny

    # Projection Kernel 
    def P_p(self, p,q, kernel='anisotropic'):

        # p = np.hstack([np.zeros(nx), p])
        # q = np.hstack([q, np.zeros(ny)])

        if kernel == 'isotropic':
            p = p / np.maximum(1, np.sqrt(p**2 + q**2) ) 
            q = q / np.maximum(1, np.sqrt(p**2 + q**2) ) 
        else: 
            p = p / np.maximum(1, np.abs(p))
            q = q / np.maximum(1, np.abs(q))

        return p, q

    # Projection onto Convex Set
    def P_c(self, x):
        return np.clip(x, 0, None)

    # Linear Operator
    def L(self, p, q):

        num_cols_p = p.shape[1]; num_rows_q = q.shape[0]
        
        p[0,:] = 0; p[-1,:] = 0
        q = np.hstack((q, np.zeros((num_rows_q, 1))))
        
        q[:,0] = 0
        p = np.vstack((p, np.zeros((1,num_cols_p))))
        
        #Produces one matrix that is m x n
        return p + q - np.roll(p, -1, axis=0) - np.roll(q, -1, axis=1)

    # Transpose of Linear Operator 
    def L_t(self, x , old_p, old_q, coeff=1):

        m,n = x.shape
        new_p = x[:-1,:] - np.roll(x[:-1,:], 1, axis=0)
        new_q = x[:,:-1] - np.roll(x[:,:-1], 1, axis=1)

        return old_p + coeff * new_p, old_q + coeff * new_q     

    # Python Implementation of FGP-TV [1] denoising / regularization (2D Case)
    # This function is based on the paper by
    # Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained Total Variation Image Denoising and Deblurring Problems"
    def fgp_tv(self, input, lambdaTV=1e2, nIter=25):  
            
        (nx, ny) = input.shape

        Px = np.zeros((nx-1, ny)); Py = np.zeros((nx, ny-1))
        Rx = Px.copy(); Ry = Py.copy()
        t = 1

        for i in range(nIter):
            
            x = self.P_c(input - lambdaTV * self.L(Px, Py))
            
            (Px, Py) = self.L_t(x, Px, Py, coeff=1 / (8 * lambdaTV))
            (Px, Py) = self.P_p(Px, Py)
            
            t_next = (1 + np.sqrt(1 + 4 * t**2)) / 2
            Rx_new = Px + ((t - 1) / t_next) * (Px - Rx)
            Ry_new = Py + ((t - 1) / t_next) * (Py - Ry)
            Rx, Ry = Px, Py
            Px, Py = Rx_new, Ry_new
            t = t_next

        return self.P_c(input - lambdaTV * self.L(Px, Py))

    def tv(self, input, kernel='isotropic'):
          
        if kernel == 'isotropic':  
            deltaX = np.diff(input,axis=0)**2
            deltaY = np.diff(input,axis=1)**2 
            tv = np.sqrt(np.sum(deltaX) + np.sum(deltaY))
        return tv