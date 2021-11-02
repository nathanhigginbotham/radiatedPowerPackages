add_library(FieldClasses FieldClasses.cxx)
target_link_libraries(FieldClasses PUBLIC BasicFunctions SignalProcessing ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})
