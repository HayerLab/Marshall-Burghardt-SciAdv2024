# Excitable Rho dynamics control cell shape and motility by sequentially activating ERM proteins and actomyosin contractility
[Seph Marshall-Burghardt](https://orcid.org/0009-0000-7634-9768), [Rodrigo Migueles-Ramirez](https://orcid.org/0000-0002-6087-1348), [Qiyao Lin](https://orcid.org/0009-0003-6436-9237), [Nada El Baba](https://orcid.org/0000-0002-2658-8426), [Rayan Saada](), [Mustakim Umar](https://orcid.org/0009-0004-6364-7431), [Kian Mavalwala](https://orcid.org/0009-0009-9631-9852), [Arnold Hayer](https://orcid.org/0000-0001-7808-8880)

Cite as "Marshall-Burghardt *et.al.*, Science Advances, 2024"

![](GraphicalAbstract.png)

- Raw data can be downloaded [here](https://doi.org/10.20383/103.01016).
- Preprint: [bioRxiv](https://doi.org/10.1101/2023.12.19.572346), to be published as Marshall-Burghardt *et.al.*, 2024 in Science Advances.

## About
This repository contains MATLAB code used for image processing and data analysis throughout our study. It is not intended as a ready-to-use resource package and it may not work without modifications for analysis of other datasets.

## Navigating the code repository 
Certain script and function packages are used multiple times throughout the paper and have their own unique folders, whereas panel-specific scripts are found in the folder 'figure panel specific scripts and instructions'. Additional instructions for analysis related to specific figure panels can be found in the README files in 'figure panel specific scripts and instructions'. 

### FRET, cell edge velocity tracking, and window analysis 
The ratiometric FRET, cell edge velocity tracking, and window analysis code uses functions written by Sean Collins, as described in:
> Yang HW., Collins SR., & Meyer, T., "Locally excitable Cdc42 signals steer cells during chemotaxis". Nat Cell Biol 18(2), 191–201 (2016). https://doi.org/10.1038/ncb3292
 
This code also uses: 
> ANN:Approximate Nearest Neighbours Version 1.1.2. Copyright (c) 1997-2010 University of Maryland and Sunil Arya and David Mount. http://www.cs.umd.edu/~mount/ANN/

### Quantitative immunofluorescence 
The quantitative immunofluorescence code uses functions written by Mingyu Chung, as described in:
> Cappell SD., Chung M., Jaimovich A., Spencer SL., Meyer T. "Irreversible APCCdh1 inactivation underlies the point of no return for cell-cycle entry". Cell 166(1), 167-180 (2016). https://doi.org/10.1016/j.cell.2016.05.077

### Violin plots 
Violin plots for data visualization were created in part with thanks to:
> Bechtold, Bastian, 2016. Violin Plots for Matlab, GitHub Project  
> https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.455984

### MATLAB toolboxes 
Running this code requires the following MATLAB toolboxes: curve fitting, image processing, statistics and machine learning, as well as [Statistical Learning Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/12333-statistical-learning-toolbox?s_tid=FX_rc2_behav) by Dahua Lin.


