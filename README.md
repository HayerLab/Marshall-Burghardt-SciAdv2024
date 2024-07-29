# Excitable Rho dynamics control cell shape and motility by sequentially activating ERM proteins and actomyosin contractility
[Seph Marshall-Burghardt](https://orcid.org/0009-0000-7634-9768), [Rodrigo Migueles-Ramirez](https://orcid.org/0000-0002-6087-1348), [Qiyao Lin](https://orcid.org/0009-0003-6436-9237), [Nada El Baba](https://orcid.org/0000-0002-2658-8426), [Rayan Saada](), [Mustakim Umar](), [Kian Mavalwala](https://orcid.org/0009-0009-9631-9852), [Arnold Hayer](https://orcid.org/0000-0001-7808-8880)

![](GraphicalAbstract.png)

- Reference to the main article: [Marshall-Burghardt *et.al.*, 2024, Science Advances]()
- Data can be found [here]().
- Link to corresponding [GitHub repository](https://github.com/HayerLab/Marshall-Burghardt-SciAdv2024).

## Understanding file names
### Default naming convention
Using the "jobs" feature of the Nikon HCS software, the microscope takes a set number of images per well and stores the images as TIFF files in folders named using the R_C_S format where R and C are the row and the column index in the 96-well plate, respectively, and S is the index of sites. Inside each site folder, each channel data gets stored into a separate channel folder and the channel name gets appended to the file name (i.e.: 2_3_1_CFP.tiff). You will see this format in most of the data, and the code has been configured for this naming format.
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
All code used in the data processing/analysis found in the paper is present here. Certain script and function packages are used multiple times throughout the paper and have their own unique folders, whereas panel-specific scripts are found in the folder 'figure panel specific scripts and instructions'.  Step-by-step instructions listing all scripts used for each panel can be found in the README files in 'figure panel specific scripts and instructions'. 

## Use of violin plots 
Throughout the paper, violinplots were sometimes created using code from 
> Bechtold, Bastian, 2016. Violin Plots for Matlab, Github Project  
> https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.455984
## How to navigate the data
The data is organized by figures in the paper. All the code necessary to preprocess, analyze the data and generate figures is inside the Code folder. 

The source images can be found under "Data" inside each date folder and then under "TIFF Stacks".
Inside each dataset folder:
* TIFF Stacks contains the source data. 
	* Crop Schemes reports on the way cells were cropped. The file key, cropping mode and parameters used are displayed. 
	  When the crop index is not specified, the image shows each site's field of view, the crop index for each cell, and the bounding boxes used to crop them. Red = Automatic cropping. Cyan = manual cropping. 
	  When the crop index is specified, the individually-inspected cell crops are reported. Cell outline and bounding box reports on the exclusion policy applied to that cell. Yellow = Cell being inspected. Red = Discarded. Magenta (automatic cropping) or cyan (manual cropping) = Other cells currently retained.
	* Cropped folder contains individually cropped cell images.
- Objects folder contains the `.mat` file with the extracted object properties.
- Metadata folder contains details on the data properties (i.e.: pixel size).

The Results folder contains aggregated statistics across many trials or conditions.
## Data structures for reproducibility
We have made our best to make sure these results are fully reproducible. In principle, it would be possible to start with the data inside of TIFF files and go through our code to process the data and generate the same results again. 
To do so, you'll need three essential things:
- `Metadata.mat` The metadata file for each dataset. This should NOT be modified.
- `Options.mat` This file is the analysis log and contains the parameters used to generate the results. It should NOT be modified. However, running the code would overwrite it, since then the way the data would have been analyzed would have changed. If you're going to attempt to reproduce the data, please make sure to:
	1. Create a copy of `Options.mat` inside of the same folder.
	2. Rename one of the versions and leave the other one with the original name. The one with the original name will be overwritten with the new analysis records.
- `ExpSettings.mat` This file contains the pointers to paths in your computer and the indexes of the datasets currently being analyzed. It is normal that this file changes constantly. For instance, you will have to update it to specify the folder you have downloaded our data to if you want to try to reproduce our results.
Please keep in mind that if you run the code to reproduce our data, the contents of the Objects, Crop Schemes, Cropped and Metadata folders will be overwritten and modified. Files that are not overwritten will not be automatically deleted.
You may need to update the paths contained in `ExpSettings` so that they point to the right folder. In this convention, the Root path is the umbrella folder that contains all other folders. In the absence of a cloud storage service like OneDrive, the Cloud path also points to the Root path.
To easily update all the paths defined in `ExpSettings`, use the `changeComputer()` function by calling ` ExpSettings = changeComputer(ExpSettings);`. If the defined Root our Cloud paths cannot be found, select them in your machine and the rest of the paths will be automatically updated. Make sure to save `ExpSettings` and overwrite the previous version, so the changes are saved.




