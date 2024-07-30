# Excitable Rho dynamics control cell shape and motility by sequentially activating ERM proteins and actomyosin contractility

## Alternative naming convention
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