# Prepare FIB-SEM Data for Segmentation

This Matlab code performs image processing and alignment tasks on *serial section images* with the goal to prepare 3D data for segmentation.

## *Overview*

### Rationale

A script-based approach allows to tailor the workflow of processing steps to the specific needs of a given image series and flexibly react to unforeseen events. The
Matlab software is well curated and its various toolboxes allow to tackle major computational challenges. 

> The goal of the given script is to provide a scaffold for a structured approach that can flexibly react to the conditions of a given image series while at the same time
maintaining a certain degree of transparency.


<details>

<summary> Examples for variability from run to run ... </summary>

- data are recorded on various FIB platforms which use different recording strategies and distinguish in their raw data structure
- a large variety of detectors are used, sometimes multiples simultaneously in a single run
- certain detector geometries lead to image gradients
- others generate noisy data that necessitate the application of denoising algorithms to facilitate the segmentation process
- image series suffer to a different degree under imperfections, such as image-to-image scan distortions, depending on the instrument and lab conditions
- materials with strong local inhomogeneity in their milling behavior show *curtaining* artifacts

</details>


### Script Structure

> Albeit being flexible, scripting caries various disadvantages. It affords a considerable knowledge of the corresponding language by the user, overall script structure can
become quickly unclear and hard to decipher, and despite possible variability there is a large number of routine processes that are used regularly. This code tries to address
some of those issues.

- The script uses *Matlab sections* (denoted by the double-percent symbol '%%' in the text) as structuring element. Each section represents *one process step* and the
  overview of currently used processes can be conveniently seen using the 'GoTo' button in the Matlab editor.
- As a result, each Matlab script represents *one particular linear workflow* of image processing tasks. The entire analysis is compactly captured in one script file.
- A list of processes is stored in a big [template file](https://github.com/StephanKraemer/FIB-SEM/blob/main/SerialDataAlignTemplate_All.m).
  A new script is built by copying the corresponding template sections. Template scripts are provided asstarting points.
- Sections themselves are built following a strict structure. They begin with user input data, followed by definition of parameters and code for the given process,
  and finish saving results into a data structure.
- A *template section* is provided with instructions for the user which allows to add new processes while maintaining the overall structure of the script.

### Ease of use

If a given workflow can be accomplished solely with preexisting processes that are provided in a template file, the user interaction is limited to adapting the variables
defined in the `UPDATE` text block which is clearly marked. It is hoped that this approach will make the script accessible and usable to novice Matlab users. Advanced
programmers can incorporate their own functions and use the script as scaffold.

### Workflow Documentation

> The general goal of this code is to streamline the scripting process and provide at the end a well-documented workflow. 

<details>

<summary> The code tries to accomplish this task in several ways... </summary>

- Intermediate image data are saved using a naming convention that reflects the list of processes applied to the given data set on the operating system level.
- Salient results of each section are saved into a Matlab `info` structure.
- A Matlab table `ProcessingList` is generated with elements documenting the various processing steps, the dimensions of the image series as they change along the way,
names of result sub structures and time needed for the execution of each step.
- Wherever useful, analysis results are plotted. In cases where plots are generated for an entire series, avi files are used and saved into the Figure folder. Series of
images are stored in the MRC format which can be read with image processing tools such as *ImageJ*. Most single plots are at this moment only shown on the
screen and the user can decide to save those as files from the plot window.

</details>



## *Documentation*

Please see the [Documentation](https://github.com/StephanKraemer/FIB-SEM/blob/main/Documentation.pdf) for a detailed description.

<details>

<summary>Table of contents </summary>

1. Overview
    1. Why Scripting, why Matlab
    2. Script structure
    3. Ease of use
    4. Workflow documentation
    5. Version 1.0
2. Script Desing
    1. Building a linear workflow
    2. Script structure of the first session
    3. Script structure of following sessions
3. Section Design
    1. Structural elements of a section
    2. Building a new process using the template section
4. Working with a Script
    1. Script templates
    2. Script navigation and execution
    3. Handling of serial section data
    4. Starting a project
    5. Stopping a session
    6. Continue a session
    7. Starting another session
5. Data Storage
    1. Serial section data
    2. Naming convention for intermediate results
    3. Processing parameters
6. Project Organization
    1. Project folders
    2. Project folder structure
7. Standard Operating Procedures
    1. Install code
    2. Set up project folders
    3. Build first session script
    4. Start first session and build info structure
    5. Parse image files on a Crossbeam FIB
    6. Parse image files on a Helios FIB
    7. Intermittently end a session
    8. Continue a session
    9. Build continuing session
8. Conventions
    1. Coordinate system
    2. Slice indexing
    3. Bounding box of an image frame
9. Key Processes
    1. Non-local means denoising
    2. Image alignment via image-to-image cross correlation
    3. Correctino of y-dependent scan distortions
10. Demo Project
11. Acknowledgements

</details>



## *Demo Project*

For illustration purposes, a complete project folder is provided representing a multi-session workflow on data recorded on a Crossbeam FIB. 
The Data belong to a research project by Yamilex Acevedo-Sanchez and Rebecca Lamason from the MIT Biology Department

> *'An obligate intracellular bacterial pathogen forms extensive and stable interkingdom contacts with the endoplasmic reticulum'*.

It can be downloaded from the Dropbox folder [Demo Project - Acevedo-Sanchez](https://www.dropbox.com/scl/fo/9ghkrjua375nmih7tn803/h?rlkey=qdk4a0m8otk0970r108gb06ot&dl=0). 
Due to space contstraints only three consecutive images containing a region of interest are uploaded. The original 16-bit images have been reduced to 8 bit.
The Matlab scripts are as applied to the entire dataset of over 1400 images.

A PowerPoint file *Summary.pptx* (in *Align* folder) contains run parameters, notes and figures produced during the analysis, illustrating the chosen degree of denoising, 
and a comparison of yz-slices before and after drift correction. For more details see chapter 'Demo Project' in the 
[Documentation](https://github.com/StephanKraemer/FIB-SEM/blob/main/Documentation.pdf).

<details>

<summary> The project contains three sessions that distinguish in the degree of processing  </summary>

- Session 1: images are denoised, drift-corrected, binned to generate isotropic voxels, and intensity-inverted to facilitate comparison with existing TEM data
- Session 2: omits denoising and applies drift-correction to original data
- Session 3: minimal processing, only scan distortions are corrected, data are not binned

</details>



## *Installation and Usage*

For detailed instructions please see chapter 'Standard Operating Procedures' in the [Documentation](https://github.com/StephanKraemer/FIB-SEM/blob/main/Documentation.pdf). 



## *System requirements*

The script was generated with ***Matlab R2022a*** on a Dell Precision Tower 7910 workstation

<details>

<summary> System configuration  </summary>

- System: 64 bit
- Memory: 128GB
- CPU: 40 cores (Intel(R) Xeon(R) CPU E5-2699 v4 @ 2.20GHz   2.20 GHz, 2 processors)
- GPU: NVIDIA Quadro K6000
- Operating System: Windows 10

</details>

> The script in its current form does *not* make use of GPU acceleration. Wherever possible `parfor` loops are used to speed up calculations.

Image data can be handled in two different ways, depending on the amount of memory available on the computer.

- In the first approach images are successively loaded from the hard drive at each processing step, and resulting images are stored back onto the hard drive after each task.
This allows to handle large datasets with a regular PC and inspect in detail the results of intermediate processing results.
- For the second approach the entire data set is loaded into memory. This case is intended for speedy evaluation.

> ***Note*** There is currently *no* memory management. If the dataset is loaded into memory, the memory has to be large enough to hold the given number
of Matlab clones, *two* representations of the data, before and after processing, plus overhead for the actual calculations.





