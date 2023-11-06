# 2023_CL062_Tao_Bhandawat

2023_CL062_Tao_Bhandawat contains the MATLAB code that accompanies the paper: "Neurons underlying aggressive actions that are shared by both males and females in Drosophila."


## Getting Started

There are 3 repositories:
The Connectomics folder contains the R and MATLAB code used to perform connectomics analysis.
The DMDCodeBase folder contains the MATLAB code used to analyze the head fixed DMD optogenetics experiments.
The WalkingChamberCodeBase contains the MATLAB code used to analyze the freely walking behaviour in the mirror chamber.

These code was written using:
MATLAB r2022a
RStudio (R v4.1.1)

You will need the following list of add ons in MATLAB (forthcoming)
1. 
2. 

You will need to install (if not already installed) the following list packages on RStudio
1. Community packages for working with flywire and hemibrain data:
	a. natverse, reticulate, neuprintr, hemibrainr, nat, natverse, nat.jrcbrains, fafbseg
2. RConnectomicsAnalysis (note that this is just located within this folder)
3. reticulate (to interface with python)
4. R.matlab (to interface with MATLAB)
5. Plotting packages
	a. pracma, lsa, pheatmap, ComplexHeatmap, RColorBrewer, gridExtra, 
6. Other useful packages
	a. progress, bit64, pracma


To get started, please download or clone the GitHub repository
```shell
$ git clone https://github.com/bhandawatlab/2023_CL062_Tao_Bhandawat.git
```

## Usage on dataset

### Connectomics
A copy of the dataset is located in [Dropbox](https://www.dropbox.com/scl/fo/u1xk2uqkyl1nrt66s8bme/h?rlkey=nwuu94l9abm1iggh8pahhsa6g&dl=0)
Please download the data and place the subfolders into the Data folder: 2023_CL062_Tao_Bhandawat/Connectomics/RConnectomicsAnalysis/Data/

In RStudio, navigate to the Connectomics folder
1. To run the analysis regarding connections from CL062, aIPg, and pC1 neurons to DNs: CL_aIPg_pC1_toDN_analysis.R
2. To run the analysis regarding connections from CL062, aIPg, and pC1 neurons to each other: CL_aIPg_pC1_connections_analysis.R
3. To run the analysis regarding connections from aSPg to CL062, aIPg, and pC1 neurons: aSP_to_CL_aIPg_pC1_analysis.R

### WalkingChamberCodeBase
A copy of the dataset is located in [Dropbox](https://www.dropbox.com/scl/fo/ncvnp7mdgky209xts1918/h?rlkey=vj6izdum2f91ry61gyyffbabw&dl=0)
Please download the data and place the subfolders into the Data folder: 2023_CL062_Tao_Bhandawat/WalkingChamberCodeBase/Data/

In MATLAB, navigate to the WalkingChamberCodeBase folder and run: RunPipeline.m

### DMDCodeBase
A copy of the dataset is located in [Dropbox](https://www.dropbox.com/scl/fo/mrb8kk4do4w3978crhwkq/h?rlkey=rqgapxeakzw9kndago72m509c&dl=0)
Please download the data and place the subfolders into the Data folder: 2023_CL062_Tao_Bhandawat/DMDCodeBase/Data/

In MATLAB, navigate to the DMDCodeBase folder and run: RunPipeline.m

Note 1: The MATLAB code creates multipage postscript analysis figure files. To convert postscript to pdf, you can use adobe or [ghostscript](https://www.ghostscript.com/)
Note 2: We will add an automatic download script for getting the data into the appropriate folders in the future.


## Other Data location

The source dataset (specific for each figure panel) can be found in [To be added](https://www.dropbox.com/scl/fo/kia94jeadfe2d2bw5sodo/h?rlkey=nf72mq1kyav3ptvkt2q3be110&dl=0)
The JAABA model can be found here [JAABA](https://www.dropbox.com/scl/fo/ypb4gwupdknzajj8qwa0m/h?rlkey=71m4rdq5j02axo7ig366ixdj5&dl=0)
The DeepLabCut models can be found here [DLC](https://www.dropbox.com/scl/fo/v6vs5h4w8ls7fjdzuz56a/h?rlkey=6x3bxfsgb7fqi0f2il4uo3a9t&dl=0)


## General Information


## Authors and Citation

Code by Liangyu Tao [lt532@drexel.edu](mailto:lt532@drexel.edu)

