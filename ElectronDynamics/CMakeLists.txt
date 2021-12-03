add_library(ElectronDynamics BaseField.cxx QTNMFields.cxx BorisSolver.cxx)
target_link_libraries(ElectronDynamics PUBLIC BasicFunctions ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})
