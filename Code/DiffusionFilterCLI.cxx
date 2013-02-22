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
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkDiscreteGaussianImageFilter.h>

// command line parsing
#include "vul_arg.h"

#include "Types.h"

DiffusionImageType::Pointer dwiImage;

bool ReadDWI( const std::string &file )
{
  typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
  FileReaderType::Pointer reader = FileReaderType::New();
  reader->SetFileName( file.c_str() );
  try
    {
    reader->Update();
    dwiImage = reader->GetOutput();
    }
  catch (itk::ExceptionObject e)
      {
      std::cout << "oops" << e << std::endl;
      return false;
      }

  return true;
}

bool WriteDWI( const std::string &file )
{
  typedef itk::ImageFileWriter<DiffusionImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( file.c_str() );
  writer->SetInput(dwiImage);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject e)
      {
      std::cout << "oops" << e << std::endl;
      return false;
      }

  return true;
}

void FilterDWI(double fwhm)
{
  DiffusionImageType::RegionType region = dwiImage->GetLargestPossibleRegion();
  DiffusionImageType::PointType origin = dwiImage->GetOrigin();
  DiffusionImageType::SpacingType spacing = dwiImage->GetSpacing();
  DiffusionImageType::DirectionType direction = dwiImage->GetDirection();

  MDImageType::Pointer tempImage;
  tempImage = MDImageType::New();
  tempImage->SetRegions( region );
  tempImage->SetOrigin( origin );
  tempImage->SetSpacing( spacing );
  tempImage->SetDirection( direction );
  tempImage->Allocate();

  typedef itk::ImageRegionIterator<DiffusionImageType> DiffusionIteratorType;
  typedef itk::ImageRegionIterator<MDImageType> MDIteratorType;

  DiffusionIteratorType dwiIt(dwiImage, dwiImage->GetLargestPossibleRegion() );
  MDIteratorType tempIt(tempImage, tempImage->GetLargestPossibleRegion() );

  for(unsigned int componentIndex = 0;
      componentIndex < dwiImage->GetVectorLength();
      ++componentIndex)
    {
    // readout into a temp image
    for ( dwiIt.GoToBegin(), tempIt.GoToBegin(); !dwiIt.IsAtEnd(); ++dwiIt, ++tempIt)
      {
      DiffusionImageType::PixelType dwiVec = dwiIt.Get();
      double value = dwiVec[componentIndex];
      tempIt.Set(value);
      }

    // filter
    typedef itk::DiscreteGaussianImageFilter< MDImageType, MDImageType >
      FilterType;
    FilterType::Pointer filter = FilterType::New();
    FilterType::ArrayType variance;
    variance[0] = fwhm*fwhm;
    variance[1] = fwhm*fwhm;
    variance[2] = fwhm*fwhm;
    filter->SetVariance(variance);
    filter->SetInput(tempImage);
    filter->Update();

    // write back filtered result
    MDIteratorType resultIt(filter->GetOutput(), filter->GetOutput()->GetLargestPossibleRegion() );
    for ( dwiIt.GoToBegin(), resultIt.GoToBegin(); !dwiIt.IsAtEnd(); ++dwiIt, ++resultIt)
      {
      DiffusionImageType::PixelType dwiVec = dwiIt.Get();
      dwiVec[componentIndex] = resultIt.Get();
      dwiIt.Set(dwiVec);
      }
    }
}

int main(int argc, char** argv)
{
  // command line args
  vul_arg<std::string> dwi_infile("-i", "Input DWI NRRD File", "");
  vul_arg<std::string> dwi_outfile("-o", "Output DWI NRRD File", "");
  vul_arg<double> filter_fwhm("-s", "Filter FWHM (mm)", 0.5);
  vul_arg_parse(argc, argv);

  if ( !ReadDWI( dwi_infile() ) ) return EXIT_FAILURE;

  FilterDWI( filter_fwhm() );

  if ( !WriteDWI( dwi_outfile() ) ) return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
