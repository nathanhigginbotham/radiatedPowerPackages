cmake_minimum_required(VERSION 3.12)
project(radiatedPower VERSION 0.1.0)

# Compiler/build settings
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Paths for FFTW ROOT wrapper
set(ENV{LIBROOTFFTW_PATH} /usr/local/lib/)
list(APPEND CMAKE_PREFIX_PATH $ENV{LIBROOTFFTW_PATH})

# Locate the ROOT package
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)
include(${ROOT_USE_FILE})
include_directories(${CMAKE_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_CXX_FLAGS})

# Add the FFTW ROOT library
find_library(ROOTFFTW_LIB RootFftwWrapper REQUIRED) 

# Add libraries
add_subdirectory(BasicFunctions)
list(APPEND EXTRA_LIBS BasicFunctions)
list(APPEND EXTRA_INCLUDES "${PROJECT_SOURCE_DIR}/BasicFunctions")

# Add any executables
add_executable(dipoleRadiationPattern dipoleRadiationPattern.cxx)
target_link_libraries(dipoleRadiationPattern BasicFunctions ${ROOT_LIBRARIES})

