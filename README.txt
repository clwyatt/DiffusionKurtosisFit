Diffusion Kurtosis Fitting and Scalar Computation
Chris Wyatt, Virginia Tech
clwyatt@vt.edu

This code consists of three programs.

Diffusion Kurtosis Fitting
--------------------------

dkifit - takes multiple DWI volumes in Nrrd format and produces
a B0 (baseline) image, DTI volume, and Diffusion Kurtosis Tensor
(DKI) volume. The tensor volumes only include the unique terms 
in the symmetric tensors.

Usage: dkifit [-i string list] [-b string] [-d string] [-k string] [-v bool] 

REQUIRED:

Optional:
  Switch Type        Help [default value] 

      -i string list Input NRRD DWI Files  []
      -b string      Output B0 File  ['B0_output.nii.gz']
      -d string      Output Diffusion Tensor File  ['DTI_output.nii.gz']
      -k string      Output Diffusion Kurtosis Tensor File  ['DKI_output.nii.gz']
      -v bool        Write progress output to stdout.  [not set]
      -h bool        Print this message

Diffusion Scalars Computation
-----------------------------
dtiscalars - computes MD and FA from the DTI output of dkifit.

Usage: dtiscalars [-i string] [-m string] [-f string] 

REQUIRED:

Optional:
  Switch Type   Help [default value] 

      -i string Input DTI File  ['']
      -m string Output Mean Diffusivity File  ['MD_output.nii.gz']
      -f string Output Fractional Anisotropy File  ['FA_output.nii.gz']
      -h bool   Print this message

Diffusion Kurtosis Computation
------------------------------
dkiscalars - compute the mean kurtosis (MK) from the DKI output of dkifit
and one of the original Nrrd inputs (to get the gradient directions).

dkiscalars [-r string] [-k string] [-d string] [-m string] [-u bool] 

REQUIRED:

Optional:
  Switch Type   Help [default value] 

      -r string Input Reference DWI NRRD File  ['']
      -k string Input KTI File  ['']
      -d string Input DTI File  ['']
      -m string Output Mean Kurtosis File  ['MK_output.nii.gz']
      -u bool   Dump individual tensor components  [not set]
      -h bool   Print this message

Converting DICOM DWI to Nrrd Format
-----------------------------------
Slicer3 can be used to convert each DWI DICOM volume to the nrrd
format with embedded gradient directions, either using the GUI
or via the command line module.

Building the Code
-----------------
Dependencies: cmake > 2.8.0 and ITK Version 3.20 

0. install cmake (www.cmake.org)
1. download and build ITK (www.itk.org)
2. set ITK_DIR to the build directory of step 1
3. make build directory, run cmake, then make

