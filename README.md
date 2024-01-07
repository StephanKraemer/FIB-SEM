# Prepare FIB-SEM Data for Segmentation

This Matlab code performs image processing and alignment tasks on *serial section images* with the goal to prepare 3D data for segmentation.

## *Overview*




## *Documentation*

Please see the 

## *Demo Project*

For illustration purposes, a complete project folder is provided representing a multi-session workflow on data recorded on a Crossbeam FIB. 
The Data belong to a research project by Yamilex Acevedo-Sanchez and Rebecca Lamason from the MIT Biology Department

> *'An obligate intracellular bacterial pathogen forms extensive and stable interkingdom contacts with the endoplasmic reticulum'*.

It can be downloaded from the Dropbox folder [Demo Project - Acevedo-Sanchez](https://www.dropbox.com/scl/fo/9ghkrjua375nmih7tn803/h?rlkey=qdk4a0m8otk0970r108gb06ot&dl=0). 
Due to space contstraints only three consecutive images containing a region of interest are uploaded. The original 16-bit images have been reduced to 8 bit.
The Matlab scripts are as applied to the entire dataset of over 1400 images.

The project contains three *sessions* that distinguish in the degree of processing

- Session 1: images are denoised, drift-corrected, binned to generate isotropic voxels, and intensity-inverted to facilitate comparison with existing TEM data
- Session 2: omits denoising and applies drift-correction to original data
- Session 3: minimal processing, only scan distortions are corrected, data are not binned

A PowerPoint file *Summary.pptx* (in *Align* folder) contains run parameters, notes and figures produced during the analysis, illustrating the chosen degree of denoising, 
and a comparison of yz-slices before and after drift correction. For more details see the *Documention*.


## *Installation and Usage*

For detailed instructions please see chapter 'Standard Operating Procedures' in the *Documentation* PDF. 


## *System requirements*

The script was generated with ***Matlab R2022a*** on a Dell Precision Tower 7910 workstation

- System: 64 bit
- Memory: 128GB
- CPU: 40 cores (Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz   2.20 GHz, 2 processors)
- GPU: NVIDIA Quadro K6000
- Operating System: Windows 10

The script in its current form does *not* make use of GPU acceleration. Wherever possible `parfor` loops are used to speed up calculations.

Image data can be handled in two different ways, depending on the amount of memory available on the computer.

- In the first approach images are successively loaded from the hard drive at each processing step, and resulting images are stored back onto the hard drive after each task.
- For the second approach the entire data set is loaded into memory. This case is intended for speedy evaluation.

> ***Note*** There is currently *no* memory management. If the dataset is loaded into memory, the memory has to be large enough to hold the given number
of Matlab clones, *two* representations of the data, before and after processing, plus overhead for the actual calculations.





