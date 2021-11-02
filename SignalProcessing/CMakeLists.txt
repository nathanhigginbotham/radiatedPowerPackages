add_library(SignalProcessing NoiseFunc.cxx Spectrogram.cxx)
target_link_libraries(SignalProcessing PUBLIC ${ROOT_LIBRARIES} ${ROOTFFTW_LIB} BasicFunctions FieldClasses)
