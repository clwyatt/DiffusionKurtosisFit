/*****************************************************************************
Copyright (c) 2012, Bioimaging Systems Lab, Virginia Tech
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of Virgina Tech nor the names of its contributors may
   be used to endorse or promote products derived from this software without
   specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*******************************************************************************/
#include <iostream>
#include <string>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <NrrdIO.h>
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

// command line parsing
#include "vul_arg.h"

#include "Types.h"

// lookup.h defines index maps from symmetric tensors
// to optimization free variable x
// constructed using make_lookup.py
// 6 should be subtracted from the index returned
// since the DTI components are not present
// TODO: clean this up, it is a mess
#include "lookup.h"

TensorImageType::Pointer dtImage;
TensorImageType::Pointer ktImage;
MDImageType::Pointer mdImage;
std::vector<DiffusionEncodingDirection> encodings;

void ReadDiffusionTensorImage(std::string filename)
{
  typedef itk::ImageFileReader<TensorImageType> FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName( filename.c_str() );
  try
    {
    reader->Update();
    dtImage = reader->GetOutput();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write Diffusion Tensor Image: " << e << std::endl;
    }
}

void ReadKurtosisTensorImage(std::string filename)
{
  typedef itk::ImageFileReader<TensorImageType> FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName( filename.c_str() );
  try
    {
    reader->Update();
    ktImage = reader->GetOutput();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write Diffusion Tensor Image: " << e << std::endl;
    }
}

void ReadEncodings(std::string filename)
{

  Nrrd *nin;
  nin = nrrdNew();

  if( nrrdLoad(nin, filename.c_str(), NULL) )
    {
    // TODO cleanup first
    return;
    }

  DiffusionEncodingDirection encoding;
  int numberKeys = nrrdKeyValueSize(nin);
  for (int keyIndex = 0; keyIndex < numberKeys; ++keyIndex)
    {
    char *key, *value;
    nrrdKeyValueIndex(nin, &key, &value, keyIndex);
    std::istringstream iss(value);
    if(!strncmp(key, "DWMRI_b-value", 13))
      {
      iss >> encoding.bvalue;
      }
    if(!strncmp(key, "DWMRI_gradient_", 15))
      {
      iss >> encoding.gx >> encoding.gy >> encoding.gz;
      encodings.push_back(encoding);
      }
    free(key); free(value);
    key = NULL; value = NULL;
    }

  nrrdNuke(nin);
}

void ComputeMeanKurtosis()
{
  TensorImageType::RegionType region = ktImage->GetLargestPossibleRegion();
  TensorImageType::PointType origin = ktImage->GetOrigin();
  TensorImageType::SpacingType spacing = ktImage->GetSpacing();
  TensorImageType::DirectionType direction = ktImage->GetDirection();

  mdImage = MDImageType::New();
  mdImage->SetRegions( region );
  mdImage->SetOrigin( origin );
  mdImage->SetSpacing( spacing );
  mdImage->SetDirection( direction );
  mdImage->Allocate();

  typedef itk::ImageRegionIterator<TensorImageType> TensorIteratorType;
  typedef itk::ImageRegionIterator<MDImageType> MDIteratorType;

  TensorIteratorType ktIt(ktImage, ktImage->GetLargestPossibleRegion() );
  TensorIteratorType dtIt(dtImage, dtImage->GetLargestPossibleRegion() );
  MDIteratorType mdIt(mdImage, mdImage->GetLargestPossibleRegion() );
  for ( ktIt.GoToBegin(), dtIt.GoToBegin(), mdIt.GoToBegin(); !ktIt.IsAtEnd(); ++ktIt, ++mdIt, ++dtIt)
      {
      TensorImageType::PixelType dtVec = dtIt.Get();
      TensorImageType::PixelType ktVec = ktIt.Get();
      double meankurt = 0;
      unsigned int numNonZeroEncodings = 0;
      for(unsigned int e = 0; e < encodings.size(); ++e)
	{
	if(encodings[e].isZero()) continue;
	numNonZeroEncodings += 1;
	double g[3];
	g[0] = encodings[e].gx;
	g[1] = encodings[e].gy;
	g[2] = encodings[e].gz;
	double dapp = 0;
	for(unsigned int i = 0; i < 3; ++i)
	  for(unsigned int j = 0; j < 3; ++j)
	  {
	  dapp += g[i]*g[j]*dtVec[D[i][j]];
	  }
	double kapp = 0;
	for(unsigned int i = 0; i < 3; ++i)
	  for(unsigned int j = 0; j < 3; ++j)
	    for(unsigned int k = 0; k < 3; ++k)
	      for(unsigned int l = 0; l < 3; ++l)
		{
		kapp += (g[i]*g[j]*g[k]*g[l]*ktVec[K[i][j][k][l]-6]);
		}
	meankurt += kapp;
	}
      double value = meankurt/static_cast<double>(numNonZeroEncodings);
      if (value > 500) value = 500;
      if (value < 0) value = 0;

      mdIt.Set(value);
      }
}

void WriteMeanKurtosisImage(std::string filename)
{
  typedef itk::ImageFileWriter<MDImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(mdImage);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write MD Image: " << e << std::endl;
    }
}

void WriteKurtosisComponents(std::string basename)
{

  TensorImageType::RegionType region = ktImage->GetLargestPossibleRegion();
  TensorImageType::PointType origin = ktImage->GetOrigin();
  TensorImageType::SpacingType spacing = ktImage->GetSpacing();
  TensorImageType::DirectionType direction = ktImage->GetDirection();

  MDImageType::Pointer tempImage;
  tempImage = MDImageType::New();
  tempImage->SetRegions( region );
  tempImage->SetOrigin( origin );
  tempImage->SetSpacing( spacing );
  tempImage->SetDirection( direction );
  tempImage->Allocate();

  typedef itk::ImageRegionIterator<TensorImageType> TensorIteratorType;
  typedef itk::ImageRegionIterator<MDImageType> MDIteratorType;

  TensorIteratorType ktIt(ktImage, ktImage->GetLargestPossibleRegion() );
  MDIteratorType tempIt(tempImage, tempImage->GetLargestPossibleRegion() );

  for(unsigned int componentIndex = 0; componentIndex < 15; ++componentIndex)
    {
    for ( ktIt.GoToBegin(), tempIt.GoToBegin(); !ktIt.IsAtEnd(); ++ktIt, ++tempIt)
      {
      TensorImageType::PixelType ktVec = ktIt.Get();
      double value = ktVec[componentIndex];

      // if(value > 10) value = 10;
      // if(value < -10) value = -10;

      tempIt.Set(value);
      }

    typedef itk::ImageFileWriter<MDImageType> FileWriterType;
    FileWriterType::Pointer writer = FileWriterType::New();
    std::ostringstream oss;
    oss << basename.c_str() << "_" << componentIndex << ".nii.gz";
    writer->SetFileName( oss.str().c_str() );
    writer->SetInput(tempImage);
    try
      {
      writer->Update();
      }
    catch (itk::ExceptionObject e)
      {
      std::cout << "Error: Could not Write MD Image: " << e << std::endl;
      }
    }
}

int main(int argc, char** argv)
{
  // command line args
  vul_arg<std::string> reffile(0, "Input Reference DWI NRRD File", "");
  vul_arg<std::string> dki_infile(0, "Input KTI File", "");
  vul_arg<std::string> dti_infile(0, "Input DTI File", "");
  vul_arg<std::string> md_outfile("-m", "Output Mean Kurtosis File", "MK_output.nii.gz");
  vul_arg<bool> dump_components("-u", "Dump individual tensor components", false);
  vul_arg_parse(argc, argv);

  ReadEncodings( reffile() );
  ReadDiffusionTensorImage( dti_infile() );
  ReadKurtosisTensorImage( dki_infile() );
  ComputeMeanKurtosis();
  WriteMeanKurtosisImage( md_outfile() );

  if( dump_components() )
    {
    WriteKurtosisComponents(std::string("DKI_component"));
    }

  return EXIT_SUCCESS;
}
