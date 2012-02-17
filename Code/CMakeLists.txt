
ADD_EXECUTABLE( dkifit
Optimizer.cxx
DiffusionKurtosisFittingApp.cxx
DiffusionKurtosisFittingCLI.cxx)

TARGET_LINK_LIBRARIES ( dkifit vul ${ITK_LIBRARIES})
