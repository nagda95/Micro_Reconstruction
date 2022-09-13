# MICROSTRUCTURE RECONSTRUCTION
**An open source framework for 3D Characterization and Reconstruction of Heterogeneous Battery Electrode Microstructure**\
This repository contains the framework for statistical-equivalent reconstruction of battery electrode microstructure from the tomography image data. The high cost involved in 3D microscopy imaging can be averted in an economical way using this stochastic reconstruction framework to derive large numbers of virtual 3D microstructures from limited morphological information of 2D cross-sectional images.

![Micro_Reconstruction](https://github.com/nagda95/Micro_Reconstruction/blob/main/Segmented_Data/recons_github.png)

This framework has been developed to support lithium-ion battery modelling, for high fidelity microstructure-resolved physics-based 3D battery model to understand the interplay between microstructure and battery performance. 

**The framework has been submitted for publication in an international peer-reviewed journal article.**

## Installation/Execution/Implementation
This framework has been tested and is compatible with MATLAB R2020b and all the above versions.  

Copy the folder ‘root’ from the Github repository in your computer and add it to the MATLAB path, including all the files and subfolders.  

The framework also requires running a R script under MATLAB, so install the R package and add the R executive file location to the MATLAB path.  

Run the MATLAB file, reconstruction_framework.m for the complete execution of the framework of stochastic reconstruction of battery electrode microstructure.  

**Tomography data:**
ETH Zürich hosts tomographic datasets for NMC-based Porous Electrodes and they are available open source for download at http://dx.doi.org/10.5905/ethz-iis-1 from ETH Zürich library, a member oft the DataCite initiative (www.datacite.org).

**Machine Learning segmentation:**
We used Stardist for Machine Learning based segmentation of active material particles in the 2D tomography images. The Stardist model can be implemented by downloading the code from the following Github repository:
https://github.com/stardist/stardist.git


## How to cite
If you use this framework, or use some or parts of the algorithms contained within the framework, please quote them accordingly:
*	If you are characterizing particle shape morphology using spherical harmonics, then please quote: Feinauer, J., Stochastic 3D modeling and simulation of electrode materials in lithium-ion batteries. 2019, Universität Ulm.
* For microstructure regeneration, please quote: Mollon, G. and J. Zhao, 3D generation of realistic granular samples based on random fields theory and Fourier shape descriptors. Computer Methods in Applied Mechanics and Engineering, 2014. 279: p. 46-65. https://doi.org/10.1016/j.cma.2014.06.022     
* If you use the tomography data from ETH Zürich, then please also quote: Ebner, M., et al., X‐ray tomography of porous, transition metal oxide based lithium ion battery electrodes. Advanced Energy Materials, 2013. 3(7): p. 845-850.

If you have suggestions, feedbacks, or want to report a bug, please contact the author at nagda@kth.se or drop a message in the discussion section of this repository.
