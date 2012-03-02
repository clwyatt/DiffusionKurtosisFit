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

  void WriteB0Image(std::string filename);
  void WriteDiffusionTensorImage(std::string filename);
  void WriteKurtosisTensorImage(std::string filename);
  void WriteResidualImage(std::string filename);

  void SetVerbosity(bool verbose)
    {
      m_verbose = verbose;
    }
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
  MDImageType::Pointer m_ResidualImage;

  bool m_verbose;
};

