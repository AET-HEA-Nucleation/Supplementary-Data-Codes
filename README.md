# Supplementary Data Codes 

**Crystal nucleation and growth in high-entropy alloys revealed by atomic electron tomography**

Yakun Yuan<sup>1,2</sup>, Saman Moniri<sup>1*</sup>, Yao Yang<sup>1*</sup>, Jihan Zhou<sup>1*</sup>, Andrew Yuan<sup>1*</sup>, Dennis S. Kim<sup>1*</sup>, Yongsoo Yang<sup>1*</sup>, Chenyang Li<sup>3</sup>, Wei Chen<sup>3</sup>, Peter Ercius<sup>4</sup>, Jianwei Miao<sup>1†</sup>*

*<sup>1</sup>Department of Physics and Astronomy and California NanoSystems Institute, University of California, Los Angeles, CA 90095, USA.*    
*<sup>2</sup>Future Material Innovation Center, School of Materials Science and Engineering and Zhangjiang Institute for Advanced Study, Shanghai Jiao Tong University, Shanghai 200240, China.*     
*<sup>3</sup>Department of Materials Design and Innovation, University at Buffalo, The State University of New York, Buffalo, NY, 14260, USA.*     
*<sup>4</sup>National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, Berkeley, CA 94720, USA.*    
*†Correspondence and requests for materials should be addressed to J.M. (j.miao@ucla.edu).*     

## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Repositary Contents](#repositary-contents)

# Overview

Nucleation in high-entropy alloys (HEAs) involves the three-dimensional atomic structural evolution and chemical ordering of multi-principal elements, presenting significant challenges to modern experimental and theoretical approaches. Using atomic electron tomography, we determine the 3D atomic structure and chemical composition of 7,662 HEA nuclei and 498 medium-entropy alloy nuclei at different stages of nucleation. We observe that the nuclei exhibit local structural order that gradually decreases from cores to boundaries and correlates with local chemical order. As the nuclei grow, the local structural order improves. At the late stage of nucleation, most nuclei coalesce without misorientation, while a small fraction form coherent twin boundaries. To explain our experimental observations, we develop an equation called gradient nucleation pathways, where the nucleation energy barrier gradually increases through multiple evolving intermediate states. These findings advance our understanding of nucleation and growth in HEAs and offer valuable insights into the underlying nucleation mechanisms. The experimental data and source codes for the 3D image reconstruction and post analysis are provided here.

# System Requirements

## Hardware Requirements

We recommend a computer with 16G DRAM, standard i7 4-core CPU, and a GPU to run most data analysis source codes. But for the 3D reconstruction of the experimental data with RESIRE, atomic tracing and refinement, we recommend a computer with large memory (256G DRAM, 16-core CPU and 1 GPU).

## Software Requirements

### OS Requirements

This package has been tested on the following Operating System:

Linux: CentOS 6 2.6.32    
Windows: Windows 10 18368.778    
Mac OSX: We have not tested it on a Mac yet, but it should in principle work.     

### Matlab Version Requirements

This package has been tested with `Matlab` R2023b. All the codes have to run in their own folders. We recommend the use of `Matlab` version R2023b or higher to test the data and source codes.

# Repositary Contents

### 1. Experiment Data

Folder: [1_Measured_data](./1_Measured_data)

This folder contains experimental images after denoising and alignment as well as their corresponding tilt angles for the 25 high entropy nanoparticles (named HEA1 - HEA25).

### 2. The REal Space Iterative REconstruction (RESIRE) Package

Folder: [2_RESIRE_package](./2_RESIRE_package)

Run the codes `Main_RESIRE_HEA.m` to perform the 3D reconstruction of the HEA nanoparticles.
### 3. Reconstructed 3D Volume

Folder: [3_Final_reconstruction_volume](./3_Final_reconstruction_volume)

This folder contains the 3D reconstructed volumes of the 25 HEA nanoparticles.

### 4. Atom_tracing_and_classification

Folder: [4_Atom_tracing_and_classification](./4_Atom_tracing_and_classification)

Run the codes `Main_atom_tracing_and_classification_HEA.m` to trace the potential atomic positions from the reconstructed 3D volumes, separate non-atoms from the potential atoms using the K-mean clustering method, and perform local classification of elemental species into type 1-3. By carefully comparing the individual atomic positions in the potential atomic models with the 3D reconstructions, a small fraction of unidentified or misidentified atoms were manually corrected, producing the 3D atomic models of the 25 HEA nanoparticles.

### 5. Atomic Position Refinement

Folder: [5_Position_refinement](./5_Position_refinement)

Run the codes `Main_position_refinement_HEA_nanoparticle.m` to refine the 3D atomic models of the 25 HEA nanoparticles.

### 6. Experimental Atomic Model

Folder: [6_Final_coordinates](./6_Final_coordinates)

This folder contains the final 3D atomic models of the 25 HEA nanoparticles.

### 7. Post Data Analysis

Folder: [7_Data_analysis](./7_Data_analysis)

Run the codes `Main_1_pdf_and_voronoi_HEA_nanoparticle.m` to calculate and plot the pair distribution functions and Voronoi analysis of the HEA nanoparticles.

Run the codes `Main_2_BOO_evolution_in_nuclei_HEA_nanoparticle.m` to compute and plot the bond orientation order parameters and identified nuclei as well as their BOO evolutions found in the HEA nanoparticles.

Run the codes `Main_3_CSRO_BOO_HEA_nanoparticle.m` to quantify the chemical short range order parameters for the HEA nanoparticles and their correlation with crystallinity.

Run the codes `Main_4_nuclei_coalescence_HEA_nanoparticle.m` to analyze the correlation between misorientation and distance of nuclei pairs in the HEA nanoparticles, and identify different types of twins formed during nuclei coalescence.# Supplementary-Data-Codes
 
