


# 2023_CL062_Tao_Bhandawat

2023_CL062_Tao_Bhandawat contains the code that accompanies the paper: 

Tao L, Ayembem D, Barranca VJ, Bhandawat V (2024) Neurons underlying aggression-like actions that are shared by both males and females in *Drosophila*. The Journal of Neuroscience.doi: https://doi.org/10.1523/JNEUROSCI.0142-24.2024

## Getting Started
**There are 3 folders:**
1. The DMDCodeBase folder contains the MATLAB code used to analyze the head fixed DMD optogenetics experiments.
2. The WalkingChamberCodeBase folder contains the MATLAB code used to analyze the freely walking behaviour in the mirror chamber.
3. The Connectomics folder contains the R and MATLAB code used to perform connectomics analysis.


**These code was written using:**\
MATLAB r2022a\
RStudio (R v4.1.1)

**You will need to install the following list of  add-ons in MATLAB**
1. Curve Fitting Toolbox
2. Image Processing Toolbox
3. Signal Processing Toolbox
4. Statistics and Machine Learning Toolbox


**You will need to install the following list of packages on RStudio**
1. Community packages for working with flywire and hemibrain data:\
	a. natverse, reticulate, neuprintr, hemibrainr, nat, natverse, nat.jrcbrains, fafbseg
2. reticulate (to interface with python)
3. R.matlab (to interface with MATLAB)
4. Plotting packages:\
   	a. pracma, lsa, pheatmap, ComplexHeatmap, RColorBrewer, gridExtra
5. Other useful packages:\
	a. progress, bit64, pracma

To get started, please download or clone the GitHub repository
```shell
$ git clone https://github.com/bhandawatlab/2023_CL062_Tao_Bhandawat.git
```

To download the data, please navigate to the **2023_CL062_Tao_Bhandawat** directory and run: ```DownloadData.m```
If this doesn't work, then please manually download the data and place it in the correct folders (see below)

### WalkingChamberCodeBase
A copy of the dataset is located in [Dropbox](https://www.dropbox.com/scl/fo/grnpfn4341rtk8lxigkoy/AHP8lHGSyFLWGSNExqhiNjw?rlkey=0maj0z4btn67hcmh9tnvol2vo&st=clm9fxwt&dl=0)\
Make sure the data is in the following folder: 2023_CL062_Tao_Bhandawat/WalkingChamberCodeBase/Data/

In MATLAB, navigate to the WalkingChamberCodeBase folder and run: ```RunPipeline.m```

### DMDCodeBase
A copy of the dataset is located in [Dropbox](https://www.dropbox.com/scl/fo/3yfi7e0bt6tk5vof2s66x/AG14bXNnMynMY_ECSjiqAFM?rlkey=zu98bvvhaj1maxy4tbedt598t&st=76j2d6ta&dl=0)\
Make sure the data is in the following folder: 2023_CL062_Tao_Bhandawat/DMDCodeBase/Data/

In MATLAB, navigate to the DMDCodeBase folder and run: ```RunPipeline.m```

### Connectomics
A copy of the dataset is located in [Dropbox](https://www.dropbox.com/scl/fo/slv34h2i3k16h6vf75ofd/AKThNByCWVBsk0D0a02_Hb0?rlkey=uq8zcrf3nhp5f6ty51klinzyv&st=wljpf6kn&dl=0)\
Make sure the data is in the following folder: 2023_CL062_Tao_Bhandawat/Connectomics/Data/

In RStudio, navigate to the Connectomics folder
1. To run the analysis regarding connections from CL062, aIPg, and pC1 neurons to DNs: ```CL_aIPg_pC1_toDN_analysis.R```
2. To run the analysis regarding connections from CL062, aIPg, and pC1 neurons to each other: ```CL_aIPg_pC1_connections_analysis.R```
3. To run the analysis regarding connections from aSPg to CL062, aIPg, and pC1 neurons: ```aSP_to_CL_aIPg_pC1_analysis.R```

In Matlab, navigate to the Connectomics folder
1. To run the PCA and SVM analysis regarding connections from CL062, aIPg, and pC1 neurons to DNs: ```Run_PCA_DN_wrapper.m```


## Other Data location

The JAABA model can be found here [JAABA](https://www.dropbox.com/scl/fo/083iqfmx6topwnqa1tfwa/AC2Uhc8d1Kh0D4putnH77SM?rlkey=5cdjfao8ukaeci9xjs2uol58a&st=pxuwlfsk&dl=0)\
The DeepLabCut models can be found here [DLC](https://www.dropbox.com/scl/fo/2wnbgbhsde420a4oyg33h/ADNc0Ks2SX-4k0VI2tj61W8?rlkey=vim53i566pynno08gvp2mumdb&st=t6z6t3ja&dl=0)


## General Information

**Note 1:** The MATLAB code creates multipage postscript analysis figure files. To convert postscript to pdf, you can use adobe or [ghostscript](https://www.ghostscript.com/)\
<br/>
**Note 2:** The figures for each of the WalkingChamberCodeBase, DMDCodeBase, and Connectomics analysis will be located in their respective Figures subfolder.

## Authors

Code by Liangyu Tao [lt532@drexel.edu](mailto:lt532@drexel.edu)

