---
title: Introduction
numbering:
  enumerator: 2.%s
---

Scanning transmission electron microscopes (STEM) can enable the precise quantification of structure and chemistry through the collection and analysis of elastic and inelastic scattering events. 
The modern electron beam can be tuned to lateral dimensions less than than an angstrom and scanned between or across atoms. 
Collecting inelastic scattering events and at every energy loss using either energy dispersive x-ray (EDX) or electron energy loss spectroscopy (EELS) has allowed routine measurement of chemistry at the atomic length-scale. 
However, the acquisition of high-resolution chemical images often necessitates doses that may surpass the tolerance of the specimen, leading to chemical maps that are either noisy or completely absent [@Hart_2017].
Challenges in this field of high resolution chemical imaging remain prevalent despite advancements in detector effiency [@McMullan_2014; @Zaluzec_2021; @Goodge_2020].
The HAADF detector collects elastically scattered electrons, offering a direct view of the atomic structure with a higher SNR, but provides limited chemical information. 
Some chemical information is inferred from HAADF projections since the elastic scattering off an element is roughly proportional to atomic number $Z^2$ [@Crewe_1970; @Hartel_1996; @LeBeau_2008; @Hovden_2012]. 

However, solving chemistry from HAADF is often difficult or impossible when mulitple chemistries are present, two elements are close in atomic number, or there are large variations in material density. 
Additionally, since the image taken on the microscope is a projection through the whole viewing area, differential density further complicates how to interpret different sections in the image and leads to false conclusions.
Traditionally, detector outputs like HAADF and EDX/EELS have been analyzed independently, overlooking the potential insights from integrating structure and chemical data from the two seperate modalities. 
Unlike correlative imaging, data-fusion merges sparse signals to enhance the accuracy of both measurements [@Hall_1997; @Lahat_2015; @Di_2016]. 
Due to its power to improve image quality, this data fusion approach has seen an explosion in popularity in the past decade in the panchromatic, infrared, hyperspectral and Lidar satellite remote sensing communities [@Sun_2023].  
Data fusion relies on crafting a model that can capture interactions between different datasets without forcing unrelated connections leading to false conclusions[@Calhoun_2016].

:::{figure} ./figs/Figure_0_comparison.png
:name: fig_overview
:width: 700px
HAADF, EDX, and Fused Multi-Modal comparison of final images collected of DyScO$_3$
:::

We walkthrough the steps for multi-modal electron microscopy as a mathematical model, a computational algorithm, and finally through code implementation to make this technique more accessible to the greater electron microscopy community.  
We have provided a sample dataset that is accompanied by a full code tutorial, and we have provided an algorithm section where every variable is explained and sequentially manipulated for full transparency. 
This approach has been tested extensively in the case on 2D projections utilizing EDX and EELS signals [@Schwartz_2022]. 
When multi-modal data fusion performs well, it significantly boosts the SNR of chemical maps by 300-500% and enables a reduction in required doses by more than an order of magnitude. 
When looking at known synthetic materials we have consistently found accurate chemical mapping and since we can determine relative concentrations within specimens we also have tested stoichiometry and found <15% error with no known priors of inelastic cross sections. 
Through this work we hope that within minutes you can go from raw data off the microscope to accurate, publishable fused multi-modal chemical maps.