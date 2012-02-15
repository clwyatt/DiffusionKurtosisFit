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
using std::cout;
using std::cerr;
using std::clog;
using std::endl;

#include <string>
#include <cassert>

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkOrientedImage.h>

// command line parsing
#include "vul_arg.h"

typedef short PixelType;
typedef itk::OrientedImage<PixelType, 3> Image3DType;
typedef itk::OrientedImage<PixelType, 2> Image2DType;

int main(int argc, char** argv)
{
  // command line args
  vul_arg<std::string> infile(0, "Input File");
  vul_arg<std::string> outfile(0, "Output File");
  vul_arg_parse(argc, argv);

  Image3DType::Pointer input;
  itk::ImageFileReader< Image3DType >::Pointer reader =
    itk::ImageFileReader< Image3DType >::New();
  reader->SetFileName(infile());

  itk::ImageFileWriter< Image3DType >::Pointer writer =
    itk::ImageFileWriter< Image3DType >::New();
  writer->SetFileName(outfile());

  writer->SetInput( reader->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Exception thrown while writing.\n";
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
