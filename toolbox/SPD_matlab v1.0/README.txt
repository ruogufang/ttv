README for Sparse Perfusion Deconvolution (SPD) Package

Ruogu Fang
Advanced Multimedia Laboratory
Department of Electrical and Computer Engineering
Cornell University
April 15, 2014


Sparse perfusion deconvolution package uses a sparse representation and maximum a posterior optimization to robustly estimate the perfusion parameters in low-dose CT perfusion. In this package, we estimate the perfusion parameter cerebral blood flow (CBF), as described in the following papers. We provide a MATLAB implementation of the KSVD-SPD and Online-SPD algorithms described in these two papers. 

Note that this implementation is an advanced version of the algorithm in the papers. Instead of use matrix inversion to optimize the second step in the iterative algorithm, we adopted a steepest descent algorithm to approach the solution with improved performance. 

CITATION:
--------------
Please cite the following papers if you use code in this SPD package.

1. Ruogu Fang, Tsuhan Chen, Pina Sanelli. Towards Robust Deconvolution of Low-Dose Perfusion CT: Sparse Perfusion Deconvolution Using Online Dictionary Learning. Medical Image Analysis, Volumn 17, Issue 4, Pages 417-428, 2013. 

2. Ruogu Fang, Tsuhan Chen, Pina Sanelli. Sparsity-Based Deconvolution of Low-Dose Perfusion CT Using Learned Dictionaries. MICCAI'12, The 15th Annual International Conference on Medical Image Computing and Computer Assisted Intervention, 2012. Lecture Notes in Computer Science Volume 7510, 2012, pp 272-280.


BIBTEX:
-----------------
@incollection{fang2012dictionary
year={2012},
isbn={978-3-642-33414-6},
booktitle={Medical Image Computing and Computer-Assisted Intervention – MICCAI 2012},
volume={7510},
series={Lecture Notes in Computer Science},
editor={Ayache, Nicholas and Delingette, Hervé and Golland, Polina and Mori, Kensaku},
title={Sparsity-Based Deconvolution of Low-Dose Perfusion CT Using Learned Dictionaries},
publisher={Springer Berlin Heidelberg},
author={Fang, Ruogu and Chen, Tsuhan and Sanelli, PinaC.},
pages={272-280}
}

@article{Fang2013417,
title = "Towards robust deconvolution of low-dose perfusion CT: Sparse perfusion deconvolution using online dictionary learning ",
journal = "Medical Image Analysis ",
volume = "17",
number = "4",
pages = "417 - 428",
year = "2013",
issn = "1361-8415",
author = "Ruogu Fang and Tsuhan Chen and Pina C. Sanelli",
}

FILES ORGANIZATION:
----------------------------------
Main_SPD.m: main file to run the demo. It tests three algorithms: cTSVD, KSVD-SPD and Online-SPD. 
cTSVD: algorithm is included in the Utilities/pct package. 
ksvdpct: KSVD-SPD algorithm which learns a dictionary using KSVD and calls ksvd_map to estimate the low-dose perfusion map.
ksvd_map: Maximum A Posterior algorithm calls ksvd_prior and ksvd_sd.
ksvd_prior: Estimate the prior image using the learned dictionary using omp solver.
spd_sd: Estimate the low-dose perfusion map image using steepest descent based on the prior image and temporal convolution model.
spd: Online SPD algorithm which calls spd_map.
spd_map: Maximum A Posterior algorithm calls spd_prior and spd_sd.
spd_prior: Estimate the prior image using the learned dictionary using lasso solver.


INSTALLATION: 
--------------------------
1.  Unzip the package. 
2.  Compile Utility packages. The binaries for Mac OS 10.9 is already included in the package. For different platforms, 
   a. Go to ompbox10 folder, follow the readme.txt file to make files in MATLAB.
   b. Go to spams-matlab folder and compile following the instructions in HOW_TO_INSTALL.txt. 
   * Note that for newer version of Mac OS (10.7+) or gcc (XCode), issues may arise when running compile.m. Please follow the notes at the end of this README file for possible solutions.
3. Run Main_SPD.m for the demo.
4. You may change the noise level and other parameters at the beginning of Main_SPD.m file. You may also use your own simulated or real low-dose image by loading different DICOM or MAT files as Vn (the low-dose CTP data [T x X x Y]). 


TOOLBOX INCLUDED:
-------------------------------
This package already includes the utility software packages downloaded from other website. Please properly cite the related papers if you use these packages.
a. metrix-mux    http://foulard.ece.cornell.edu/gaubatz/metrix_mux/
b. ompbox v10    http://www.cs.technion.ac.il/~ronrubin/software.html
c. spams-matlab v2.3  http://spams-devel.gforge.inria.fr



=========================================================================================================
NOTES for Issues at Compilation of SPAMS-MATLAB package. 
-------------------
1. Issue may occur when compiling spams_matlab package.
When compile.m is run, the error occurs:
/Applications/MATLAB_R2013a.app/extern/include/tmwtypes.h:819:9: error: unknown
      type name 'char16_t'
typedef char16_t CHAR16_T;

Reason: The problem is due to the update of OS system and gcc (or XCode), which changes the version number and system parameters, but MATLAB does not know it.

Solution:  In /Applications/MATLAB_R2013a.app/bin/mexopts.sh
On Line 150: Add -std=c++11 to CXXFLAGS in your current mexopts.sh. 
Re-run mex -setup to have MATLAB reconfigure it for you and it works.

When it is not possible to use C++11 via CXXFLAGS, either due to compiler limitation or because the file must be C only, the possible workarounds:
Add before #include "mex.h" in .c file:
#include <stdint.h>
#define char16_t UINT16_T


2. Issue: mex: link of ' "./build//mexTrainDL.mexmaci64"' failed.

Solution: Change the use_multithread to false in compile.m file


Please contact rf294@cornell.edu if you have issues or suggestions for this SPD package.
