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
#include <vector>
#include <string>

// command line parsing
#include "vul_arg.h"

#include "Types.h"
#include "DiffusionKurtosisFittingApp.h"

int main(int argc, char** argv)
{
  // command line args
  vul_arg<std::vector< std::string > > infiles("-i", "Input NRRD DWI Files");
  vul_arg<std::string> b0_outfile("-b", "Output B0 File", "B0_output.nii.gz");
  vul_arg<std::string> dti_outfile("-d", "Output Diffusion Tensor File", "DTI_output.nii.gz");
  vul_arg<std::string> dki_outfile("-k", "Output Diffusion Kurtosis Tensor File", "DKI_output.nii.gz");
  vul_arg<bool> verbose("-v", "Write progress output to stdout.", false);
  vul_arg_parse(argc, argv);

  DiffusionKurtosisFittingApp app;
  app.SetVerbosity( verbose() );

  app.ReadEncodings(infiles());
  if(!app.ReadDWI(infiles()) )
    {
    std::cout << "Error: DWI images cannot be read or are incompatible" << std::endl;
    return EXIT_FAILURE;
    }
  if(verbose()) std::cout << "Input Completed." << std::endl;

  app.ComputeB0Image();
  if(verbose()) std::cout << "B0 Computation Completed." << std::endl;

  app.ComputeDiffusionAndKurtosis();

  app.WriteB0Image(b0_outfile());
  app.WriteDiffusionTensorImage(dti_outfile());
  app.WriteKurtosisTensorImage(dki_outfile());
  app.WriteResidualImage(std::string("fit_residual.nii.gz"));
  return EXIT_SUCCESS;
}
