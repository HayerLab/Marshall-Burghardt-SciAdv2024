# Excitable Rho dynamics control cell shape and motility by sequentially activating ERM proteins and actomyosin contractility
[Seph Marshall-Burghardt](https://orcid.org/0009-0000-7634-9768), [Rodrigo Migueles-Ramirez](https://orcid.org/0000-0002-6087-1348), [Qiyao Lin](https://orcid.org/0009-0003-6436-9237), [Nada El Baba](https://orcid.org/0000-0002-2658-8426), [Rayan Saada](), [Mustakim Umar](), [Kian Mavalwala](https://orcid.org/0009-0009-9631-9852), [Arnold Hayer](https://orcid.org/0000-0001-7808-8880)

![](GraphicalAbstract.png)

- Reference to the main article: [Marshall-Burghardt *et.al.*, 2024, Science Advances]()
- Data can be found [here](https://doi.org/10.20383/103.01016).
- The latest version of our code related to this study is available on https://github.com/HayerLab/Marshall-Burghardt-SciAdv2024


## Understanding file names
### Default naming convention
Using the "jobs" feature of the Nikon HCS software, the microscope takes a set number of images per well and stores the images as TIFF files in folders named using the R_C_S format where R and C are the row and the column index in the 96-well plate, respectively, and S is the index of sites. This is followed by the channel name or abbreviation (i.e.: 2_3_1_CFP.tiff). You will see this format in most of the data, and the code has been configured for this naming format.
### Alternative naming convention
For certain figures, including Fig. 7 and Fig. S8c, we used a different naming convention.  Files acquired using the jobs feature were grouped together under a single folder and renamed. The corresponding codes for these sections use the following convention: 
Each file key is YYMMDD-FF-SSS-CC-WW-CNDTN-TT-CHNL where:
* YYMMDD: Acquisition date
* FF: File index
* SSS: Site index
* CC: Crop index (when applicable)
* WW: Well index (Ex.: B3)
* CNDTN: Well name (experimental condition nickname)
* TT: Treatment duration (in hours)
* CHNL: Channel name

## Navigating the repository 
MATLAB image analysis code used for image processing and analysis throughout our study is provided here. Certain script and function packages are used multiple times throughout the paper and have their own unique folders, whereas panel-specific scripts are found in the folder 'figure panel specific scripts and instructions'.  Additional instructions for analysis related to specific figure panels can be found in the README files in 'figure panel specific scripts and instructions'. 

## FRET, cell edge velocity tracking, and window analysis 
The ratiometric FRET, cell edge velocity tracking, and window analysis code uses functions written by Sean Collins, described in:
> Yang, HW., Collins SR., & Meyer, T., "Locally excitable Cdc42 signals steer cells during chemotaxis". Nat Cell Biol 18(2), 191–201 (2016). https://doi.org/10.1038/ncb3292
 
This code also uses: 
> ANN:Approximate Nearest Neighbours Version 1.1.2. Copyright (c) 1997-2010 University of Maryland and Sunil Arya and David Mount. http://www.cs.umd.edu/~mount/ANN/

## Quantiative immunofluorescnce 
The quantitative immunofluorescence code uses functions written by Mingyu Chung, first described in:
> Cappell SD., Chung M., Jaimovich A., Spencer SL., Meyer T. "Irreversible APCCdh1 inactivation underlies the point of no return for cell-cycle entry". Cell 166(1), 167-180 (2016). https://doi.org/10.1016/j.cell.2016.05.077

## Violin plots 
Violin plots for data visualization were created in part with thanks to:
> Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project  
> https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.455984



