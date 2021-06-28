# NCCENTROID3D

This package includes custom MATLAB scripts for 3D centroid analysis for in vivo negative contrast imaging as described in:

Wu, Jung et al. Intravital fluorescence microscopy with negative contrast. 

The scripts were last revised and tested in MATLAB R2018a, using a 3.1GHz Intel Core i7 MacBookPro with 16Gb memory. The package also includes a sample negative image z-stack and the required mask z-stacks for running a demo, as well as a copy of the DEMO output.

## Instructions for running Demo

1. Download package, place in MATLAB path. Add path.
2. Enter path .../MATLAB/NCCENTROID3D
3. Run scripts in the following order:

* VesselEst_01.m
* VesselLeak_02.m
* VesselLkClean_03.m
* BkgdEstCutoff_04.m
* CellCtSize2D_05.m
* CellCtSize3D_06.m

Each script was written to automatically draw the required input from /INPUT and /OUTPUT from the previous scripts.

## Description

NCCENTROID3D takes an z-stack of negative contrast images and locates the centroid of cells in 2D and 3D. It also calculates cell sizes and 3D distances from the cell centroids to their nearest structure of interest (vessel wall).

![*Overview of NCCENTROID3D*](https://github.com/juwellwwu/NCCENTROID3D/blob/main/README%20Img/Package%20Description%20Illustration.png?raw=true)

In addition to the negative contrast image z-stack, the codes also accept multiple image masks stacks for input. Masked areas are excluded from cell centroid analysis. The vessel mask z-stack is also used for 3D-distance-to-vessel-wall measurement. These image mask z-stacks can be created in any image processing software. Bone and vessel mask z-stacks used in this DEMO were generated in FIJI by segmenting the original second-harmonic-generation bone z-stack (blue channel in INPUT/DEMO_ImgStack_RGB) and the original FITC-dextran channel (green channel in INPUT/DEMO_ImgStack_RGB; the same data is stored in INPUT/DEMO_ImgStack_G.mat). While not used in this DEMO, the codes also allow two additional mask z-stacks ("ImgMaskStack" and "RegionMask") that allow users to block out specific 3D ("ImgMaskStack") or 2D ("RegionMask") areas of the negative contrast image z-stack from centroid anlaysis. Bone mask z-stack, vessel mask z-stack, and "ImgMaskStack" should be identical in size to the negative contrast image z-stack, while "RegionMask" should be a single image of the same height and width as the negative contrast image z-stack. Alignment of the supplied mask z-stacks and the negative contrast image z-stack is verifiable in OUTPUT/CellCtSize2D_Out/Colocal WSMap OrgImgEC, which are colour-coded based on how the code interprets each pixel for centroid analysis (yellow=cells; blue=bone; red=vessels; magenta=too-dim-to-analyze areas; green=blocked-out pixels.)
