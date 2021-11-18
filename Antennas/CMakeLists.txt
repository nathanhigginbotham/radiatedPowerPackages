add_library(Antennas IAntenna.cxx HertzianDipole.cxx HalfWaveDipole.cxx)
target_link_libraries(Antennas PUBLIC ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})
