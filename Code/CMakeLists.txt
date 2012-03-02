
ADD_EXECUTABLE( dkifit
Optimizer.cxx
DiffusionKurtosisFittingApp.cxx
DiffusionKurtosisFittingCLI.cxx
ProgressMeter.cxx)

TARGET_LINK_LIBRARIES ( dkifit vul ${ITK_LIBRARIES})

ADD_EXECUTABLE( dtiscalars
DiffusionScalarsCLI.cxx)

TARGET_LINK_LIBRARIES ( dtiscalars vul ${ITK_LIBRARIES})

ADD_EXECUTABLE( dkiscalars
KurtosisScalarsCLI.cxx)

TARGET_LINK_LIBRARIES ( dkiscalars vul ${ITK_LIBRARIES})

ADD_EXECUTABLE( dwisim
DiffusionSimulateCLI.cxx)

TARGET_LINK_LIBRARIES ( dwisim vul ${ITK_LIBRARIES})
