cmake_minimum_required(VERSION 3.18)
project(radiatedPower VERSION 0.1.0)

# Compiler/build settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Locate the ROOT package
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Add the FFTW ROOT library
find_library(ROOTFFTW_LIB RootFftwWrapper REQUIRED)
message("ROOTFFTW_LIB is set to ${ROOTFFTW_LIB}")

# Add libraries
add_subdirectory(BasicFunctions)
list(APPEND EXTRA_LIBS BasicFunctions)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/BasicFunctions")

add_subdirectory(FieldClasses)
list(APPEND EXTRA_LIBS FieldClasses)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/FieldClasses")

add_subdirectory(SignalProcessing)
list(APPEND EXTRA_LIBS SignalProcessing)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/SignalProcessing")

# Add any executables
add_executable(sineTest sineTest.cxx)
target_link_libraries(sineTest PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(dipoleRadiationPattern dipoleRadiationPattern.cxx)
target_link_libraries(dipoleRadiationPattern PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(normTest normTest.cxx)
target_link_libraries(normTest PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(spectrogramTest spectrogramTest.cxx)
target_link_libraries(spectrogramTest PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(noiseTests noiseTests.cxx)
target_link_libraries(noiseTests PUBLIC ${EXTRA_LIBS} ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

add_executable(polarPlots polarPlots.cxx)
target_link_libraries(polarPlots PUBLIC BasicFunctions ${ROOT_LIBRARIES})
