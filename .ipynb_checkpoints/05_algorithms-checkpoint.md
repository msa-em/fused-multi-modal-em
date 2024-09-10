---
title: Algorithms
numbering:
  enumerator: 6.%s
---

+++ {"part": "Algorithm"} 

### Algorithms

:::{prf:algorithm} Multi-Modal 2D Data Fusion
:label: Multi-Modal 2D Data Fusion

**Objective Function:**      
$\Psi(x) = \frac{1}{2}\left\|\sum_i x_i - b_H\right\|^2 + \lambda_{Chem}\sum_i\left(1^T x_i - b_{i}^T\log(x_i + \epsilon)\right) + \lambda_{TV}\sum_i \left\|x_i\right\|_{TV}$

**Inputs:**  
$b_H \in \mathbb{R}^{n_x \times n_y} \rightarrow$ HAADF Micrograph (Image Dimensions: $n_x \times n_y$)<br>
$b_C \in \mathbb{R}^{n_e \times n_x \times n_y} \rightarrow$ Raw Chemical Maps (Total \# of Elements: $n_e$)<br>
$N_{iter} = 30 \rightarrow$ Number of Main Cost Function Iterations<br>
$N_{iter TV} = 3 \rightarrow$ Number of Total Variation (TV) Iterations<br>
$\epsilon = 0.2 \rightarrow$ Background Noise Threshold<br>
$\lambda_{HAADF} = \frac{1}{n_e} \rightarrow$ HAADF Weight<br>
$\lambda_{Chem} = 0.08 \rightarrow$ Data Consistency Weight, Ranges from 0 to 1<br>
$\lambda_{TV} = 0.15 \rightarrow$ Regularization Weight, Ranges from 0 to 1<br>
$γ = 1.6 \rightarrow$ Incoherent Linear Imaging factor

**Output:**    
$x \in \mathbb{R}^{n_e \times n_x \times n_y} \rightarrow$ Concatenated Vector of Recovered Maps

**Main Function:**    
1. Initialize first iterate as raw elemental maps $x_i^0 = b_C$<br>
2. For $k = 1$ to $N_{iter}$:
	1. $x^k = x^{k-1} - (\lambda_{HAADF}\nabla \Psi_1(x^k) + \lambda_{Chem} \nabla \Psi_2(x^k))$
	2. For $i = 1$ to $n_e$:
		1. $x_i^k = \text{TV\_FGP\_2D}(x_i^k, \lambda_{TV}, N_{iter TV}) \text{  or  TV\_GP\_2D}(x_i^k, \lambda_{TV}, N_{iter TV})$
   3. End for
3. End for
4. Return $x^{N_{iter}}$

**Definitions:**    
$\nabla \Psi_1(x) = -\gamma\text{diag}(x^{-1})A^T(b_H - Ax^{\gamma})$<br>
$\nabla \Psi_2(x) = 1 - b \oslash (x + \epsilon)$<br> 
$\oslash$ is element wise division

:::

:::{prf:algorithm} 2D Total Variation Fast Gradient Projection Method
:label: (TV_FGP_2D)
**Inputs:**  
$b \in \mathbb{R}^{n_x \times n_y} \rightarrow$ 2D Image, $\lambda \rightarrow$ Regularization Parameter, $N_{iter TV} \rightarrow$ Number of Iterations

**Output:**    
$x \in \mathbb{R}^{n_e \times n_x \times n_y} \rightarrow$ Concatenated Vector of Recovered Maps

**TV_FGP_2D Function:**    
1.  $({r}_x, {r}_y) = ({p}_x, {p}_y) = \left({0}_{(m-1) \times n}, {0}_{m \times (n-1)}\right), \quad t_1 = 1$
2.  For $k = 1, \ldots, N_{iter TV}$:
	1. $({p}_x, {p}_y) = {P}_p \left[ ({r}_x, {r}_y) + \frac{1}{8\lambda} \mathcal{L}^T\left({P}_c\left[b - \lambda \mathcal{L}({r}_x, {r}_y)\right]\right) \right]$
    2. $t_{k+1} = \frac{1 + \sqrt{1 + 4t_k^2}}{2}$
    3. $({r}_x^{k+1}, {r}_y^{k+1}) = ({p}_x, {p}_y) + \frac{t_k - 1}{t_{k+1}}\left(({p}_x - {p}_x^{k-1}, {p}_y - {p}_y^{k-1})\right)$
3. End for
4. Return $x^* = {P}_c \left[b - \lambda \mathcal{L}({p}_x^N, {p}_y^N)\right]$
:::

:::{prf:algorithm} 2D Total Variation Gradient Projection Method
:label: (TV_GP_2D)
**Inputs:**  
$b \in \mathbb{R}^{n_x \times n_y} \rightarrow$ 2D Image, $\lambda \rightarrow$ Regularization Parameter, $N_{iter TV} \rightarrow$ Number of Iterations

**Output:**    
$x \in \mathbb{R}^{n_e \times n_x \times n_y} \rightarrow$ Concatenated Vector of Recovered Maps

**TV_GP_2D Function:**    
1.  $p_x^0 = \mathbf{0} \in \mathbb{R}^{(m-1) \times n}$, $p_y^0 = \mathbf{0} \in \mathbb{R}^{m \times (n-1)}$
2.  For $k = 1, \ldots, N_{iter TV}$:
	1. $(p_x^k, p_y^k) = P_p(p_x^{k-1}, p_y^{k-1}) + \frac{1}{8}\mathcal{L}^T(P_c(b - \lambda \mathcal{L}(p_x^{k-1}, p_y^{k-1})))$
3. End for
4. Return $x^* = P_c(b - \lambda \mathcal{L}(p_x^{N_{iter TV}}, p_y^{N_{iter TV}}))$
:::

:::{prf:algorithm} The Linear Operation 
$$\mathcal{L}: \mathbb{R}^{(m-1) \times n} \times \mathbb{R}^{m \times (n-1)} \rightarrow \mathbb{R}^{m \times n}$$

**$\mathcal{L}(p, q)_{i,j}$ Function:**    
1. For $i = 1, \ldots, n_x$
	1. For $j = 1, \ldots, n_y$
        1. If $i = 0$ or $i = n_x$
              1. $p_{i,j} = 0$
        2. End if
        3. If $j = 0$ or $j = n_y$
              1. $q_{i,j} = 0$
        4. End if
        5. $\mathcal{L}(p, q)_{i,j} = p_{i,j} + q_{i,j} - p_{i-1,j} - q_{i,j-1}$
    2. End for
2. End for
3. Return $\mathcal{L}(p, q)_{i,j}$
:::

:::{prf:algorithm} The Adjoint Linear Operation
$${L}^T: \mathbb{R}^{m \times n} \rightarrow \mathbb{R}^{(m-1) \times n} \times \mathbb{R}^{m \times (n-1)}$$
**$\mathcal{L}^T(x)$ Function:**    
1. For $i = 1, \ldots, n_x - 1$
	1. For $j = 1, \ldots, n_y - 1$
        1. $p_{i,j} = x_{i,j} - x_{i+1,j}$
        2. $q_{i,j} = x_{i,j} - x_{i,j+1}$
    2. End for
2. End for
3. Return $(p, q)_{i,j}$
:::

:::{prf:algorithm} The Orthogonal Projection Operator onto set $C$
**$P_c(x)$ Function:**    
1. Return $\max\{0, x\}$
:::

:::{prf:algorithm} The Projection onto set $P$
$$P_{\mathcal{P}}(p, q) = (r, s) \rightarrow r \in \mathbb{R}^{(m-1) \times n}, s \in \mathbb{R}^{m \times (n-1)}$$
**$P_p(p, q)$ Function:**    
1. For $i = 1, \ldots, n_x - 1$
	1. For $j = 1, \ldots, n_y - 1$
        1. $p_{i,j} = \max\{1, (p_{i,j}^2 + q_{i,j}^2)^{-\frac{1}{2}}\}p_{i,j}$
        2. $q_{i,j} = \max\{1, (p_{i,j}^2 + q_{i,j}^2)^{-\frac{1}{2}}\}q_{i,j}$
    2. End for
2. End for
3. Return $(p, q)_{i,j}$
:::