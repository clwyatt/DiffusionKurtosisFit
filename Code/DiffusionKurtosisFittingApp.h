#include <vector>
#include <string>

#include "Types.h"

class DiffusionKurtosisFittingApp
{
public:

  DiffusionKurtosisFittingApp();

  bool ReadEncodings(std::vector< std::string > files);

  bool ReadDWI(std::vector< std::string > files);

  bool CheckCompatibilityDWI();

  void ComputeB0Image();

  void ComputeDiffusionAndKurtosis();

  void ComputeMeanDiffusion();

  void WriteB0Image(std::string filename);
  void WriteMDImage(std::string filename);

  void PrintInfo();

private:

  void CollapseDWI();
  void ComputeNonZeroEncodings();

  unsigned int m_NumberZeroEncodings, m_NumberNonZeroEncodings;
  std::vector<DiffusionImageType::Pointer> m_dwiImages;
  std::vector<DiffusionEncodingDirection> m_dwiEncodings;
  DiffusionImageType::Pointer m_FullEncodingImage;
  TensorImageType::Pointer m_DiffusionTensorImage;
  TensorImageType::Pointer m_KurtosisTensorImage;

  B0ImageType::Pointer m_B0Image;
  MDImageType::Pointer m_MDImage;
};

