#include <vector>
#include <string>

#include "Types.h"

class DiffusionKurtosisFittingApp
{
public:

  DiffusionKurtosisFittingApp();

  bool ReadEncodings(std::vector< std::string > files);

  bool ReadDWI(std::vector< std::string > files);

  void ComputeDiffusionAndKurtosis();

  void WriteKurtosisTensorImage(std::string filename);
  void WriteDiffusionTensorImage(std::string filename);

private:

  bool CheckCompatibilityDWI();
  void CollapseDWI();
  void AllocateResult();

  // input
  DiffusionImageType::RegionType m_region;
  DiffusionImageType::PointType m_origin;
  DiffusionImageType::SpacingType m_spacing;
  DiffusionImageType::DirectionType m_direction;
  std::vector<DiffusionImageType::Pointer> m_dwiImages;
  std::vector<DiffusionEncodingDirection> m_dwiEncodings;

  // processing
  unsigned int m_NumberVoxels;
  double *dwiData;

  // output
  TensorImageType::Pointer m_DiffusionTensorImage;
  TensorImageType::Pointer m_KurtosisTensorImage;

};
