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
#include <itkImageRegionIterator.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

// command line parsing
#include "vul_arg.h"

#include "Types.h"

TensorImageType::Pointer dtImage;
MDImageType::Pointer mdImage;
MDImageType::Pointer faImage;

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


void ComputeMeanDiffusion()
{
  TensorImageType::RegionType region = dtImage->GetLargestPossibleRegion();
  TensorImageType::RegionType::IndexType index = region.GetIndex();
  TensorImageType::RegionType::SizeType size = region.GetSize();
  TensorImageType::PointType origin = dtImage->GetOrigin();
  TensorImageType::SpacingType spacing = dtImage->GetSpacing();
  TensorImageType::DirectionType direction = dtImage->GetDirection();

  mdImage = MDImageType::New();
  mdImage->SetRegions( region );
  mdImage->SetOrigin( origin );
  mdImage->SetSpacing( spacing );
  mdImage->SetDirection( direction );
  mdImage->Allocate();

  faImage = MDImageType::New();
  faImage->SetRegions( region );
  faImage->SetOrigin( origin );
  faImage->SetSpacing( spacing );
  faImage->SetDirection( direction );
  faImage->Allocate();

  typedef itk::ImageRegionIterator<TensorImageType> TensorIteratorType;
  typedef itk::ImageRegionIterator<MDImageType> MDIteratorType;

  TensorIteratorType dtIt(dtImage, dtImage->GetLargestPossibleRegion() );
  MDIteratorType mdIt(mdImage, mdImage->GetLargestPossibleRegion() );
  MDIteratorType faIt(faImage, faImage->GetLargestPossibleRegion() );
  for ( dtIt.GoToBegin(), mdIt.GoToBegin(), faIt.GoToBegin(); !dtIt.IsAtEnd(); ++dtIt, ++mdIt, ++faIt)
      {
      TensorImageType::PixelType dtVec = dtIt.Get();
      vnl_matrix<double> DT(3,3);
      DT(0,0) = dtVec[0];
      DT(0,1) = dtVec[1];
      DT(0,2) = dtVec[2];
      DT(1,0) = dtVec[1];
      DT(1,1) = dtVec[3];
      DT(1,2) = dtVec[4];
      DT(2,0) = dtVec[2];
      DT(2,1) = dtVec[4];
      DT(2,2) = dtVec[5];

      vnl_symmetric_eigensystem<double> E(DT);
      double l1 = E.get_eigenvalue(0);
      double l2 = E.get_eigenvalue(1);
      double l3 = E.get_eigenvalue(2);

      double md = (l1 + l2 + l3)/3.0;
      mdIt.Set(md);

      double fanum = sqrt(3*((l1-md)*(l1-md) + (l2-md)*(l2-md) + (l3-md)*(l3-md)));
      double faden = sqrt(2*(l1*l1 + l2*l2 + l3*l3));
      double fa = fanum/faden;
      faIt.Set(fa);
      }
}

void WriteMeanDiffusionImage(std::string filename)
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

void WriteFractionalAnisotropyImage(std::string filename)
{
  typedef itk::ImageFileWriter<MDImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(faImage);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write FA Image: " << e << std::endl;
    }
}


int main(int argc, char** argv)
{
  // command line args
  vul_arg<std::string> infile("-i", "Input DTI File", "");
  vul_arg<std::string> md_outfile("-m", "Output Mean Diffusivity File", "MD_output.nii.gz");
  vul_arg<std::string> fa_outfile("-f", "Output Fractional Anisotropy File", "FA_output.nii.gz");
  vul_arg_parse(argc, argv);

  ReadDiffusionTensorImage( infile() );
  ComputeMeanDiffusion();
  WriteMeanDiffusionImage( md_outfile() );
  WriteFractionalAnisotropyImage( fa_outfile() );

  return EXIT_SUCCESS;
}
