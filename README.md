Diffusion Kurtosis Fitting and Scalar Computation
==================================================
by Chris Wyatt, Virginia Tech (clwyatt@vt.edu)

This code consists of three programs.

Diffusion Kurtosis Fitting
--------------------------

dkifit - takes multiple DWI volumes in Nrrd format and produces
a DTI volume and a Diffusion Kurtosis Tensor (DKI) volume.
The tensor volumes only include the unique terms in the symmetric tensors
in row-major order. Uses the same approach as described in:
Tabesh, A., Jensen, J. H., Ardekani, B. A., & Helpern, J. A. (2011).
Estimation of tensors and tensor-derived measures in diffusional kurtosis imaging.
Magnetic resonance in medicine, 65(3), 823–36. doi:10.1002/mrm.22655
except it uses an interior-point method for the constrained optimization.

	Usage: dkifit [ string list] [-d string] [-k string]

	REQUIRED:
	string list Input NRRD DWI Files  []

	Optional:
	Switch Type        Help [default value]

	-d string      Output Diffusion Tensor File  ['DTI_output.nii.gz']
	-k string      Output Diffusion Kurtosis Tensor File  ['DKI_output.nii.gz']
	-h bool        Print this message

Diffusion Scalars Computation
-----------------------------
dtiscalars - computes MD and FA from the DTI output of dkifit.

	Usage: dtiscalars [ string] [-m string] [-f string]

	REQUIRED:
	string Input DTI File  ['']

	Optional:
	Switch Type   Help [default value]

	-m string Output Mean Diffusivity File  ['MD_output.nii.gz']
	-f string Output Fractional Anisotropy File  ['FA_output.nii.gz']
	-h bool   Print this message

Diffusion Kurtosis Computation
------------------------------
dkiscalars - compute the mean kurtosis (MK) from the DKI output of dkifit
and one of the original Nrrd inputs (to get the gradient directions).

	Usage: dkiscalars [ string] [ string] [ string] [-m string] [-u bool]

	REQUIRED:
	string Input Reference DWI NRRD File  ['']
	string Input KTI File  ['']
	string Input DTI File  ['']

	Optional:
	Switch Type   Help [default value]

	-m string Output Mean Kurtosis File  ['MK_output.nii.gz']
	-u bool   Dump individual tensor components  [not set]
	-h bool   Print this message

Converting DICOM DWI to Nrrd Format
-----------------------------------
Slicer 3.6 or 4.x can be used to convert each DWI DICOM volume to the nrrd
format with embedded gradient directions, either using the GUI
or via the command line module.

Building the Code
-----------------
Dependencies: cmake >= 2.8.5 and ITK Version 4.3.1
(and an ITK  supported toolchain)

0. install cmake (www.cmake.org)
1. download and build ITK (www.itk.org)
2. set ITK_DIR to the build directory of step 1
3. make build directory, run cmake, then make
4. make test runs a test on a small volume

To take advantage of multiple cores modify the compiler
flags to enable openmp for your system (e.g. -fopenmp with g++).
This runs each voxel fit in a separate thread. Enabling loop-unrolling
will get you another few ms per voxel decrease.

Example
-------
The Testing directory also contains a larger subvolume consisting
of two files: test_medium_dataset1.nrrd and test_medium_dataset2.nrrd.

To estimate the tensors run:
	dkifit test_medium_dataset1.nrrd and test_medium_dataset2.nrrd
This will produce two files: DTI_ouput.nii.gz and DKI_output.nii.gz
The number of bad fits (if any) is reported.

To compute the diffusion scalars MD and FA run:
	dtiscalars DTI_output.nii.gz
This will produce the files FA_output.nii.gz and MD_output.nii.gz

To compute the mean kurtosis run:
	dkiscalars test_small_dataset1.nrrd DKI_output.nii.gz DTI_output.nii.gz
This will produce the file MK_ouput.nii.gz

The time to fit the tensors is approximately 300 ms per voxel per cpu
with optimized code generation (cmake BUILD_TYPE=Release).
