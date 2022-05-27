add_library(Antennas IAntenna.cxx HertzianDipole.cxx HalfWaveDipole.cxx PatchAntenna.cxx)
target_link_libraries(Antennas PUBLIC ${ROOT_LIBRARIES})
