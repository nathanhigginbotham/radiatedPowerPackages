add_library(ElectronDynamics BaseField.cxx QTNMFields.cxx)
target_link_libraries(ElectronDynamics PUBLIC ${ROOT_LIBRARIES} ${Boost_MATH_LIBRARY})
