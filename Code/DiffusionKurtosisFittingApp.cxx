#include <iostream>
#include <sstream>
#include <cassert>

#include <itkImageSeriesReader.h>
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMetaDataObject.h>
#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>

#include "DiffusionKurtosisFittingApp.h"
#include "Optimizer.h"

#include <vnl/algo/vnl_symmetric_eigensystem.h>

DiffusionKurtosisFittingApp::DiffusionKurtosisFittingApp()
{
  m_NumberZeroEncodings = 0;
  m_NumberNonZeroEncodings = 0;
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
	if(encoding.isZero() ) m_NumberZeroEncodings += 1;
	else m_NumberNonZeroEncodings += 1;
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
  DiffusionImageType::RegionType::IndexType index = region.GetIndex();
  DiffusionImageType::RegionType::SizeType size = region.GetSize();
  DiffusionImageType::PointType origin = m_dwiImages[0]->GetOrigin();
  DiffusionImageType::SpacingType spacing = m_dwiImages[0]->GetSpacing();
  DiffusionImageType::DirectionType direction = m_dwiImages[0]->GetDirection();

  m_FullEncodingImage = DiffusionImageType::New();
  m_FullEncodingImage->SetRegions( region );
  m_FullEncodingImage->SetOrigin( origin );
  m_FullEncodingImage->SetSpacing( spacing );
  m_FullEncodingImage->SetDirection( direction );
  m_FullEncodingImage->SetVectorLength( m_NumberNonZeroEncodings + m_NumberZeroEncodings );
  m_FullEncodingImage->Allocate();

  typedef itk::ImageRegionIterator<DiffusionImageType> IteratorType;
  IteratorType newIt(m_FullEncodingImage, region);
  unsigned int encodingIndex = 0;
  for(unsigned int i = 0; i < m_dwiImages.size(); ++i)
    {
    IteratorType nextIt(m_dwiImages[i], m_dwiImages[i]->GetLargestPossibleRegion());
    nextIt.GoToBegin();
    unsigned int numberEncodingsThisImage = nextIt.Get().GetSize();
    for ( newIt.GoToBegin(); !newIt.IsAtEnd(); ++newIt, ++nextIt)
      {
      DiffusionImageType::PixelType newVec = newIt.Get();
      DiffusionImageType::PixelType nextVec = nextIt.Get();
      for(unsigned int j = 0; j < numberEncodingsThisImage; ++j)
	{
	if(nextVec[j] == 0) newVec[encodingIndex +j] = 1; // clip to
							  // prevent
							  // log(0) later
	else newVec[encodingIndex +j] = nextVec[j];
	}
      newIt.Set(newVec);
      }
    encodingIndex += numberEncodingsThisImage;
    }

  // release original data
  m_dwiImages.clear();

}

void DiffusionKurtosisFittingApp::ComputeB0Image()
{
  DiffusionImageType::RegionType region = m_FullEncodingImage->GetLargestPossibleRegion();
  DiffusionImageType::RegionType::IndexType index = region.GetIndex();
  DiffusionImageType::RegionType::SizeType size = region.GetSize();
  DiffusionImageType::PointType origin = m_FullEncodingImage->GetOrigin();
  DiffusionImageType::SpacingType spacing = m_FullEncodingImage->GetSpacing();
  DiffusionImageType::DirectionType direction = m_FullEncodingImage->GetDirection();

  m_B0Image = B0ImageType::New();
  m_B0Image->SetRegions( region );
  m_B0Image->SetOrigin( origin );
  m_B0Image->SetSpacing( spacing );
  m_B0Image->SetDirection( direction );
  m_B0Image->Allocate();

  typedef itk::ImageRegionIterator<DiffusionImageType> DWIIteratorType;
  typedef itk::ImageRegionIterator<B0ImageType> B0IteratorType;

  DWIIteratorType dwiIt(m_FullEncodingImage, region);
  B0IteratorType b0It(m_B0Image, m_B0Image->GetLargestPossibleRegion() );
  for ( dwiIt.GoToBegin(), b0It.GoToBegin(); !dwiIt.IsAtEnd(); ++dwiIt, ++b0It)
      {
      DiffusionImageType::PixelType dwiVec = dwiIt.Get();
      InternalPixelType average = 0;
      for(unsigned int j = 0; j < dwiVec.GetSize(); ++j)
	{
	if(m_dwiEncodings[j].isZero()) average += dwiVec[j]/m_NumberZeroEncodings;
	}
      if(average < 1) average = 1.; // clip to prevent problems with
				    // tensor computations later
      b0It.Set(average);
      }
}

void DiffusionKurtosisFittingApp::ComputeDiffusionAndKurtosis()
{
  DiffusionImageType::RegionType region = m_FullEncodingImage->GetLargestPossibleRegion();
  DiffusionImageType::RegionType::IndexType index = region.GetIndex();
  DiffusionImageType::RegionType::SizeType size = region.GetSize();
  DiffusionImageType::PointType origin = m_FullEncodingImage->GetOrigin();
  DiffusionImageType::SpacingType spacing = m_FullEncodingImage->GetSpacing();
  DiffusionImageType::DirectionType direction = m_FullEncodingImage->GetDirection();

  m_DiffusionTensorImage = TensorImageType::New();
  m_DiffusionTensorImage->SetRegions( region );
  m_DiffusionTensorImage->SetOrigin( origin );
  m_DiffusionTensorImage->SetSpacing( spacing );
  m_DiffusionTensorImage->SetDirection( direction );
  m_DiffusionTensorImage->SetVectorLength( 6 );
  m_DiffusionTensorImage->Allocate();

  m_KurtosisTensorImage = TensorImageType::New();
  m_KurtosisTensorImage->SetRegions( region );
  m_KurtosisTensorImage->SetOrigin( origin );
  m_KurtosisTensorImage->SetSpacing( spacing );
  m_KurtosisTensorImage->SetDirection( direction );
  m_KurtosisTensorImage->SetVectorLength( 15 );
  m_KurtosisTensorImage->Allocate();

  typedef itk::ImageRegionIterator<DiffusionImageType> DiffusionIteratorType;
  typedef itk::ImageRegionIterator<TensorImageType> TensorIteratorType;
  typedef itk::ImageRegionIterator<B0ImageType> B0ImageIteratorType;

  DiffusionIteratorType dwiIt(m_FullEncodingImage, m_FullEncodingImage->GetLargestPossibleRegion() );
  B0ImageIteratorType b0It(m_B0Image, m_B0Image->GetLargestPossibleRegion() );
  TensorIteratorType dtIt( m_DiffusionTensorImage,  m_DiffusionTensorImage->GetLargestPossibleRegion() );
  TensorIteratorType ktIt( m_KurtosisTensorImage,  m_KurtosisTensorImage->GetLargestPossibleRegion() );

  unsigned int max_iterations = size[0]*size[1]*size[2];
  unsigned int iterations = 1;
  for(dwiIt.GoToBegin(), b0It.GoToBegin(), dtIt.GoToBegin(), ktIt.GoToBegin();
      !dwiIt.IsAtEnd();
      ++dwiIt, ++b0It, ++dtIt, ++ktIt)
    {
    // pull out vectors for this voxel and convert to vnl type
    DiffusionImageType::PixelType dwiVec = dwiIt.Get();
    vnl_vector<double> vnl_dwi(dwiVec.GetSize());
    for(unsigned int i = 0; i < dwiVec.GetSize(); ++i) vnl_dwi[i] = dwiVec[i];

    InternalPixelType b0 = b0It.Get();

    // construct cost function
    optimizer_cost_function cost(&vnl_dwi, b0, &m_dwiEncodings);

    // initialize x
    vnl_vector<double> x(21);
    optimizer_init_x(x);

    // find optimum
    vnl_levenberg_marquardt minimizer(cost);
    if(!minimizer.minimize_without_gradient(x))
      {
      std::cout << "minimizer failed." << std::endl;
      minimizer.diagnose_outcome();
      }

    // fill Diffusion tensor after shifting negative eigenvalues
    TensorImageType::PixelType dtVec = dtIt.Get();
    for(unsigned int i = 0; i < 6; ++i)
      {
      dtVec[i] = x[i];
      }

    // shift negative eigenvalues
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

    if(E.get_eigenvalue(0) < 0)
      {
      dtVec[0] += -E.get_eigenvalue(0);
      dtVec[3] += -E.get_eigenvalue(0);
      dtVec[5] += -E.get_eigenvalue(0);
      }
    dtIt.Set(dtVec);

    // fill Kurtosis tensor
    TensorImageType::PixelType ktVec = ktIt.Get();
    for(unsigned int i = 6; i < 21; ++i)
      {
      ktVec[i-6] = x[i];
      }
    ktIt.Set(ktVec);

    //std::cout << iterations << " / " << max_iterations << std::endl;
    iterations += 1;
    }
}

void DiffusionKurtosisFittingApp::ComputeMeanDiffusion()
{
  DiffusionImageType::RegionType region = m_FullEncodingImage->GetLargestPossibleRegion();
  DiffusionImageType::RegionType::IndexType index = region.GetIndex();
  DiffusionImageType::RegionType::SizeType size = region.GetSize();
  DiffusionImageType::PointType origin = m_FullEncodingImage->GetOrigin();
  DiffusionImageType::SpacingType spacing = m_FullEncodingImage->GetSpacing();
  DiffusionImageType::DirectionType direction = m_FullEncodingImage->GetDirection();

  m_MDImage = MDImageType::New();
  m_MDImage->SetRegions( region );
  m_MDImage->SetOrigin( origin );
  m_MDImage->SetSpacing( spacing );
  m_MDImage->SetDirection( direction );
  m_MDImage->Allocate();

  typedef itk::ImageRegionIterator<TensorImageType> DWIIteratorType;
  typedef itk::ImageRegionIterator<MDImageType> MDIteratorType;

  DWIIteratorType dtIt(m_DiffusionTensorImage, m_DiffusionTensorImage->GetLargestPossibleRegion() );
  MDIteratorType mdIt(m_MDImage, m_MDImage->GetLargestPossibleRegion() );
  for ( dtIt.GoToBegin(), mdIt.GoToBegin(); !dtIt.IsAtEnd(); ++dtIt, ++mdIt)
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

      double meandiff = (E.get_eigenvalue(0) + E.get_eigenvalue(1) + E.get_eigenvalue(2))/3.0;
      mdIt.Set(meandiff);
      }
}

void DiffusionKurtosisFittingApp::WriteB0Image(std::string filename)
{
  typedef itk::ImageFileWriter<B0ImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(m_B0Image);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write B0 Image: " << e << std::endl;
    }
}

void DiffusionKurtosisFittingApp::WriteMDImage(std::string filename)
{
  typedef itk::ImageFileWriter<MDImageType> FileWriterType;
  FileWriterType::Pointer writer = FileWriterType::New();
  writer->SetFileName( filename.c_str() );
  writer->SetInput(m_MDImage);
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject e)
    {
    std::cout << "Error: Could not Write MD Image: " << e << std::endl;
    }
}

void DiffusionKurtosisFittingApp::PrintInfo()
{
  std::cout << "Number of Encodings: " << m_dwiEncodings.size() << std::endl;
  std::cout << m_NumberZeroEncodings << " are Zero " << std::endl;
  std::cout << m_NumberNonZeroEncodings << " are Non-Zero " << std::endl;

  std::cout << m_FullEncodingImage << std::endl;
}
