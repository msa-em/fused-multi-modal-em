---
title: Mathematical Overview
numbering:
  enumerator: 3.%s
---

:::{figure} ./figs/Figure_1_math.png
:name: fig_math
:width: 700px
Overview of Fused Multi-Modal pipeline on an atomic-res DyScO$_3$ STEM/X-EDS dataset
:::

### Cost function

Fused multi-modal electron microscopy solves an optimization problem to recover micro- and nanoscale material chemistry. 
Solutions to the optimization problem are sought that encourage sparsity in the gradient domain while correlating high SNR HAADF projections with elementally sensitive EELS and/or X-EDS maps. 
These solutions ensure reduced spatial variation and accurate elemental mapping. 
The overall optimization function is framed as an inverse problem taking the following form:

$$
\label{math_eq_1}
\hat{x} = \arg\min_{x\geq 0} \left( \Psi_1(x) + \lambda_1 \Psi_2(x) + \lambda_2 \text{TV}(x) \right),
$$

where $\hat{x}$ is the final reconstruction, $\Psi_1$, $\Psi_2$, $\text{TV}$ are each optimization functions, and $\lambda$ are their respective weights.

Each term is made up of standard optimization functions which is best illustrated by the full expression:

$$
\label{math_eq_2}
\arg \min_{x_i \geq 0} \frac{1}{2} \left\| b_H - \sum_i (Z_i x_i)^\gamma \right\|_2^2 
+ \lambda_1 \sum_i (1^T x_i - b_i^T \log(x_i + \epsilon))
+ \lambda_2 \sum_i \| x_i \|_{TV} 
$$

Breaking down each variable, $b_H$ is the measured HAADF, $b_i$ and $x_i$ are the measured and reconstructed chemical maps for element $i$, $\epsilon$ is a number close to 0 that can account for the background and prevents $\log(0)$, $\log$ is applied element-wise to its arguments, superscript $T$ denotes vector transpose, and $\mathbf{1}$ denotes the vector of $n_x \times n_y$ ones, where $n_x \times n_y$ is the image size. 
When implementing an algorithm to solve this problem, we concatenate the multi-element spectral variables $(x_i, b_i)$ as a single vector: $x, b \in \mathbb{R}^{n_x \times n_y \times n_i}$, where $n_i$ denotes the total number of reconstructed elements.

Equation [](#math_eq_2) can be broken down into three key terms. 
The first term is our assumption of a forward model where HAADF projections can be described as a linear combination of elemental distributions raised to the power $\gamma \in [1.4, 2]$ [@Hartel_1996]. 
Incoherent linear imaging for high-angle elastic scattering scales with atomic number Z raised to $\gamma$. 
The second term ensures that any recovered signals maintain data fidelity with the initial inputs (raw measurements).  
A maximum negative log-likelihood for elemental maps dominated by low-count Poisson statistics is employed to enforce this data fidelity [@Di_2017]. 
The final term is common channel-wise total variation ($TV$) regularization [@Rudin_1992]. 
$TV$ reduces noise by prioritizing the preservation of sharp features. 
Each of these three terms has a weight attributed that must be balanced to ensure convergence and accurate elemental recovery.

By descending along with the negative gradient directions for the first two terms and subsequently evaluating the isotropic TV proximal operator to denoise the chemical maps we can minimize the cost function [@Beck_2009]. 
The first two term gradients are:
$$
\nabla_x \Psi_1(x) = -\gamma \text{diag} \left( (x^\gamma - 1) \right) A^T \left( b_H - Ax^\gamma \right)
$$

$$
\nabla_x \Psi_2(x) = 1 - b \oslash (x + \epsilon)
$$

where $\oslash$ denotes point-wise division. 
Here, the first term in the cost function, relating the elastic and inelastic modalities, has been equivalently re-written as ${\Psi_1 = \frac{1}{2} \left\| b_H - Ax^\gamma \right\|^2}$, where ${A \in \mathbb{R}^{n_x \cdot n_y \times n_{x} \cdot n_{y} \cdot n_i}}$ expresses the summation of all elements as matrix--vector multiplication. 
Evaluation for the TV proximal operator is in itself another iterative algorithm. 
In addition, we impose a non-negativity constraint since negative concentrations are unrealistic.
We initialize the first iterate with the measured data ${(x^0_i = b_i)}$, an ideal starting point as it is a local minima for $\Psi_2$.

$\frac{1}{L}$ is the reciprocal of the Lipschitz constant $L$ and serves as a theoretical upper limit for the step size, ensuring possible convergence.
Through the principle of Lipschitz continuity, the step size related to the gradient of the model term $\nabla \Psi_1$ was determined to be no greater than $\frac{1}{(\Vert A \Vert _1 \Vert A \Vert _{\infty})} = \frac{1}{n_i}$.
The gradient associated with the Poisson negative log-likelihood, $\Psi_2$, is not Lipschitz continuous.
This means that it is impossible to select its descent parameter ahead of time [@Dupe_2009].

### Total Variation

:::{figure} ./figs/Figure_2_TV_Comparison.png
:name: fig_tv
:width: 700px
Comparison of GP and FGP convergence on a Fused Multi-Modal dataset showing how FGP converges must faster than GP
:::

TV minimization [@Rudin_1992], regularizes images while preserving sharp features, denoising, and deblurring. 
Since TV efficiently retains edges during denoising, it is an ideal candidate for electron microscope projection processing due to the low SNR environment often encountered when trying to image materials. 
TV minimization involves solving a complex optimization problem where the norm actively penalizes the roughness of candidate solutions, thus ensuring smoothness everywhere except at sharp edges. 
A regularization parameter must be selected to balance fata fidelity with the observed data to the noise inherent in the image.

Traditionally, the Gradient Projection (GP) method has been used to solve for TV minimization, but recent advancements in optimization have introduced a faster convergence method, FGP.
First introduced by [@Beck_2009], FGP boasts computational efficiency with only minimal additions to standard Gradient Projection steps.
Just like GP, FGP can handle both isotropic and anistropic TV functions making it an ideal candidate for various image processing scenarios.
In terms of convergence rates, FGP exhibits a $O\left(\frac{1}{k^2}\right)$ rate of convergence instead of $O\left(\frac{1}{k}\right)$ found in GP.
By implementing FGP we have seen at least 10x improvement in convergence speed in our code allowing for hundreds of iterations to be run in under 30 seconds on a Mac laptop.
This implementation of FGP allows for parameter fine tuning to be done in minites allowing for Fused MM results from raw data in a single short post-processing session.
