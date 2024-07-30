# Excitable Rho dynamics control cell shape and motility by sequentially activating ERM proteins and actomyosin contractility

## Understanding file names

### Default naming convention
Using the "jobs" feature of the Nikon NIS software, the microscope takes a set number of images per well and stores the images as TIFF files in folders named using the R_C_S format where R and C are the row and the column index in the 96-well plate, respectively, and S is the index of sites. This is followed by the channel name or abbreviation (i.e.: 2_3_1_CFP.tiff). You will see this format in most of the data, and the code has been configured for this naming format.

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

## Fig. S8.
### Effect of Cpd31 on nocodazole-induced Rho activity (Fig. S8 C)
open preprocessing and single cell analysis
1. `Setup.m` 
2. `Preprocessing.m`
3. `FRET.m`

return to this folder
4. `PlotNocCPD31onRhoActivity`