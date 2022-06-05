add_library(FieldClasses FieldClasses.cxx)
target_link_libraries(FieldClasses PUBLIC BasicFunctions ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})
