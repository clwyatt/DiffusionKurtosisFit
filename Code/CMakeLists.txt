
ADD_EXECUTABLE( dkifit DiffusionKurtosisFittingCLI.cxx )
TARGET_LINK_LIBRARIES ( dkifit 
vul 
ITKAlgorithms 
ITKBasicFilters 
ITKCommon
ITKIO)
