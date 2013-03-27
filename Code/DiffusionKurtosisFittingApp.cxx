#include <iostream>
#include <sstream>
#include <cassert>

#include <itkImageSeriesReader.h>
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMetaDataObject.h>
#include <NrrdIO.h>
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>

#include "DiffusionKurtosisFittingApp.h"

#include "Optimizer.h"

#include <vnl/algo/vnl_symmetric_eigensystem.h>

DiffusionKurtosisFittingApp::DiffusionKurtosisFittingApp()
{
  m_NumberVoxels = 0;
}

bool DiffusionKurtosisFittingApp::ReadEncodings(std::vector< std::string > files)
{

  std::vector< std::string >::iterator fit;
  for(fit = files.begin(); fit != files.end(); ++fit)
    {
    Nrrd *nin;
    nin = nrrdNew();

    if( nrrdLoad(nin, fit->c_str(), NULL) )
      {
      // TODO cleanup first
      return false;
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
	m_dwiEncodings.push_back(encoding);
	}
      free(key); free(value);
      key = NULL; value = NULL;
      }

    nrrdNuke(nin);
    }

  return true;
}

bool DiffusionKurtosisFittingApp::ReadDWI(std::vector< std::string > files)
{

  std::vector< std::string >::iterator fit;
  for(fit = files.begin(); fit != files.end(); ++fit)
    {
    typedef itk::ImageFileReader<DiffusionImageType> FileReaderType;
    FileReaderType::Pointer reader = FileReaderType::New();
    reader->SetFileName( fit->c_str() );
    try
      {
      reader->Update();
      m_dwiImages.push_back( reader->GetOutput() );
      }
    catch (itk::ExceptionObject e)
      {
      std::cout << "oops" << e << std::endl;
      return false;
      }
    }

  // output image geometry information is taken from the first dwi image
  m_region = m_dwiImages[0]->GetLargestPossibleRegion();
  m_origin = m_dwiImages[0]->GetOrigin();
  m_spacing = m_dwiImages[0]->GetSpacing();
  m_direction = m_dwiImages[0]->GetDirection();

  if( !CheckCompatibilityDWI() ) return false;

  CollapseDWI();

  return true;
}

bool DiffusionKurtosisFittingApp::CheckCompatibilityDWI()
{
  DiffusionImageType::RegionType region0 = m_dwiImages[0]->GetLargestPossibleRegion();
  DiffusionImageType::RegionType::IndexType index0 = region0.GetIndex();
  DiffusionImageType::RegionType::SizeType size0 = region0.GetSize();
  DiffusionImageType::PointType origin0 = m_dwiImages[0]->GetOrigin();
  DiffusionImageType::SpacingType spacing0 = m_dwiImages[0]->GetSpacing();
  bool compatible = true;
  for(unsigned int i = 1; i < m_dwiImages.size(); ++i)
    {
    DiffusionImageType::RegionType region = m_dwiImages[i]->GetLargestPossibleRegion();
    DiffusionImageType::RegionType::IndexType index = region.GetIndex();
    DiffusionImageType::RegionType::SizeType size = region.GetSize();
    DiffusionImageType::PointType origin = m_dwiImages[i]->GetOrigin();
    DiffusionImageType::SpacingType spacing = m_dwiImages[i]->GetSpacing();

    // check Index
    compatible = compatible && (index0 == index);

    // check Size
    compatible = compatible && (size0 == size);

    // check Spacing
    compatible = compatible && (spacing0 == spacing);

    // check Origin
    compatible = compatible && (origin0 == origin);
    }

  return compatible;
}

void DiffusionKurtosisFittingApp::CollapseDWI()
{
  assert(m_dwiImages.size() > 0);

  DiffusionImageType::RegionType region = m_dwiImages[0]->GetLargestPossibleRegion();
  DiffusionImageType::RegionType::SizeType size = region.GetSize();

  m_NumberVoxels = size[0]*size[1]*size[2];
  unsigned int numberEncodings = m_dwiEncodings.size();

  dwiData = new double[m_NumberVoxels*numberEncodings];

  typedef itk::ImageRegionIterator<DiffusionImageType> IteratorType;
  unsigned int encodingIndex = 0;
  for(unsigned int i = 0; i < m_dwiImages.size(); ++i)
    {
      IteratorType nextIt(m_dwiImages[i], m_dwiImages[i]->GetLargestPossibleRegion());
      nextIt.GoToBegin();
      unsigned int voxel = 0;
      unsigned int numberEncodingsThisImage = nextIt.Get().GetSize();
      for (; !nextIt.IsAtEnd(); ++nextIt)
	{
	  DiffusionImageType::PixelType nextVec = nextIt.Get();
	  for(unsigned int j = 0; j < numberEncodingsThisImage; ++j)
	    {
	      // clip to prevent log(0) later
	      double data = (nextVec[j] == 0) ? 1. : nextVec[j];
	      dwiData[voxel*numberEncodings + encodingIndex +j] = data;
	    }
	  ++voxel;
	}
      encodingIndex += numberEncodingsThisImage;
    }

  // release original data
  m_dwiImages.clear();
}

void DiffusionKurtosisFittingApp::AllocateResult()
{
  m_DiffusionTensorImage = TensorImageType::New();
  m_DiffusionTensorImage->SetRegions( m_region );
  m_DiffusionTensorImage->SetOrigin( m_origin );
  m_DiffusionTensorImage->SetSpacing( m_spacing );
  m_DiffusionTensorImage->SetDirection( m_direction );
  m_DiffusionTensorImage->SetVectorLength( 6 );
  m_DiffusionTensorImage->Allocate();

  m_KurtosisTensorImage = TensorImageType::New();
  m_KurtosisTensorImage->SetRegions( m_region );
  m_KurtosisTensorImage->SetOrigin( m_origin );
  m_KurtosisTensorImage->SetSpacing( m_spacing );
  m_KurtosisTensorImage->SetDirection( m_direction );
  m_KurtosisTensorImage->SetVectorLength( 15 );
  m_KurtosisTensorImage->Allocate();
}

void DiffusionKurtosisFittingApp::ComputeDiffusionAndKurtosis()
{
  unsigned int numberEncodings = m_dwiEncodings.size();

  // allocate temp result space
  double * result = new double[21*m_NumberVoxels];

  unsigned int numberBadFits = 0;

  #pragma omp parallel shared(numberBadFits)
  {
    Optimizer opt(m_dwiEncodings);

    vnl_vector_fixed<double, 21> X;

    #pragma omp for
    for(unsigned int voxel = 0; voxel < m_NumberVoxels; ++voxel)
      {
	opt.SetDWI(&dwiData[voxel*numberEncodings]);

	// initial condition is spherical diffusion, small kurtosis
	// note some care is needed here to ensure that any constraint
	// is negative, but not too close to zero
	for(unsigned int i = 0; i < 21; i++) X[i] = 0;
	X[0] = 1; X[3] = 1; X[5] = 1;
	X[6] = 1e-5; X[16] = 1e-5; X[20] = 1e-5;

	// note: uses previous result as initial condition
	bool badfit;
	opt.solve(X, badfit);

	if(badfit)
	  {
	    #pragma omp critical
	    numberBadFits += 1;
	  }

	for(unsigned int i = 0; i < 21; ++i) result[voxel*21 + i] = X[i];
       }
  }

	    std::cout << numberBadFits << " bad fits ("
		      << static_cast<double>(numberBadFits)/static_cast<double>(m_NumberVoxels)
		      << " percent)" << std::endl;

  // free raw dwi data
  delete [] dwiData;

  // allocate and fill in result images
  AllocateResult();

  typedef itk::ImageRegionIterator<TensorImageType> TensorIteratorType;
  TensorIteratorType dtIt( m_DiffusionTensorImage,  m_DiffusionTensorImage->GetLargestPossibleRegion() );
  TensorIteratorType ktIt( m_KurtosisTensorImage,  m_KurtosisTensorImage->GetLargestPossibleRegion() );
  unsigned int voxel = 0;
  for(dtIt.GoToBegin(), ktIt.GoToBegin();
      !(dtIt.IsAtEnd() || ktIt.IsAtEnd());
      ++dtIt, ++ktIt)
    {
      vnl_vector_fixed<double, 21> X;
      for(unsigned int i = 0; i < 21; ++i) X[i] = result[voxel*21 + i];

      // fill Diffusion tensor
      TensorImageType::PixelType dtVec = dtIt.Get();
      for(unsigned int i = 0; i < 6; ++i)
	{
	  dtVec[i] = X[i];
	}

      double meanDiff = (dtVec[0]+dtVec[3]+dtVec[5])/3;

      // fill Kurtosis tensor
      TensorImageType::PixelType ktVec = ktIt.Get();
      for(unsigned int i = 6; i < 21; ++i)
	{
	  ktVec[i-6] = X[i]/(meanDiff*meanDiff);
	}
      ktIt.Set(ktVec);

      ++voxel;
    }

  // free temp result space
  delete [] result;
}


void DiffusionKurtosisFittingApp::WriteDiffusionTensorImage(std::string filename)
{
  typedef itk::ImageFileWriter<TensorImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(m_DiffusionTensorImage);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write Diffusion Tensor Image: " << e << std::endl;
    }
}

void DiffusionKurtosisFittingApp::WriteKurtosisTensorImage(std::string filename)
{
  typedef itk::ImageFileWriter<TensorImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(m_KurtosisTensorImage);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write Kurtosis Tensor Image: " << e << std::endl;
    }
}
