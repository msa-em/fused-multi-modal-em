---
title: Fused Multi-Modal Electron Microscopy
subject: Tutorial
subtitle: A beginner's guide
---

+++ {"part": "abstract"} 
Tomographic mapping of nano-scale chemistry via electron microscopy has challenged researchers for decades due to the fundamental tradeoff between electron dose, resolution, and reconstruction accuracy. Inelastic spectroscopic signals (electron energy loss (EELS) or energy-dispersive x-ray (EDX) spectroscopy) that yield elemental signatures are rare, resulting in long dwell times to produce resolvable structure. This process often destroys the sample and distorts the chemical map during acquisition.  Lower dose scans can be taken to combat this issue, but they often do not have enough signal to produce meaningful results. 

Fused multi-modal electron microscopy offers nano- and atomic-resolution high signal-to-noise-ratio (SNR) recovery of material chemistry by merging low-dose elastically scattered high-angle annular dark-field (HAADF) signals with the inelastic signals ([Schwartz 2022](https://doi.org/10.1038/s41524-021-00692-5)) ([Schwartz 2024](https://doi.org/10.1038/s41467-024-47558-0)). This fusion results in an order-of-magnitude dose reduction and can be easily performed on local computer hardware in minutes. This article bridges the gap between theory and practice by walking through each step of the fused multi-modal computational pipeline and pointing out best practices when running the code. By the end of this tutorial, anyone with inelastic and elastic 2D projection data should be able to reproducibly  reconstruct chemical maps using fused multi-modal electron microscopy.

:::{figure} ./figs/Figure_abstract.png
:name: fig_overview
:width: 700px
Fused Multi-Modal atomic resolution image improvement on  DyScO$_3$ 
:::

+++ {"part": "acknowledgements"} 

### Acknowledgements

We acknowledge support from Dow Chemical Company. 

Work at the Molecular Foundry was supported by the Office of Basic Energy Sciences, of the U.S. Department of Energy under Contract no. DE-AC02-05CH11231.

This research used resources of the Argonne Leadership Computing Facility, a U.S. Department of Energy (DOE) Office of Science user facility at Argonne National Laboratory and is based on research supported by the U.S. DOE Office of Science-Advanced Scientific Computing Research Program, under Contract No. DE-AC02-06CH11357.

This research used resources of the Oak Ridge Leadership Computing Facility at the Oak Ridge National Laboratory, which is supported by the Office of Science of the U.S. Department of Energy under Contract No. DE-AC05-00OR22725.

+++ {"part": "competing interests"} 

### Competing Interests

The authors declare that they have no known competing interests.
