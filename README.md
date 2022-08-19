# BS-IDT
This is the github repository for the paper [**Bond-selective intensity diffraction tomography**](https://arxiv.org/abs/2203.13630). This repository includes the codes and links to example data for implementing the BS-IDT method.

**Abstract**

Recovering molecular information remains a grand challenge in the widely used holographic and computational imaging technologies. To address this challenge, we developed a computational mid-infrared photothermal microscope, termed Bond-selective Intensity Diffraction Tomography (BS-IDT). Based on a low-cost brightfield microscope with an add-on pulsed light source, BS-IDT recovers both infrared spectra and bond-selective 3D refractive index maps from intensity-only measurements. High-fidelity infrared fingerprint spectra extraction is validated. Volumetric chemical imaging of biological cells is demonstrated at a speed of ~20 seconds per volume, with a lateral and axial resolution of ~350 nm and ~1.1 µm, respectively. BS-IDT’s application potential is investigated by chemically quantifying lipids stored in cancer cells and volumetric chemical imaging on Caenorhabditis elegans with a large field of view (~100 µm x 100 µm).

**Requirements**

This code was implemented in Matlab 2022a. When operating in other matlab versions, please be advised that some aspects of this code may not function as expected. 

**Running the Codes**

There are two main scripts for this system. The first, Main.m, runs the BS-IDT reconstruction code. The second, Main_DPC.m, runs the Bond-Selective Differential Phase Contrast (BS-DPC) reconstruction code also used in this work. The BS-IDT code uses the data as-is, while the BS-DPC code synthesizes half circle annular illuminations from the BS-IDT system and applies the reconstruction code from the work "Quantitative differential phase contrast imaging in an LED array microscope" by Tian and Waller to reconstruct 2D phase maps of the object. Additional information for running these scripts can be found within the code. 

1. Download this github repository to your computer.
2. Download the example data from our [Google Drive](https://drive.google.com/drive/folders/15Jb9bYLUzktRw0cosDe07ZOzHYORQkt4?usp=sharing)
3. Move the downloaded data to the 'data' folder within the same path as the BS-IDT code. 
4. Download any additional matlab packages required for operating the code. 

** Notes **

This code utilizes the angle self-calibration code originally developed by [Laura Waller's Lab](https://github.com/Waller-Lab/Angle_SelfCalibration) with an additional nonlinear fitting algorithm developed in the work "High-speed *in vitro* intensity diffraction tomography" By Li and Matlock et al. This prior work was used as the baseline code for BS-IDT and can be found [here](https://github.com/bu-cisl/IDT-using-Annular-Illumination).


