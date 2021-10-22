cmake_minimum_required(VERSION 3.18)
project(radiatedPower VERSION 0.1.0)

# Compiler/build settings
set(CMAKE_CXX_STANDARD 14)
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

# Add any executables
add_executable(dipoleRadiationPattern dipoleRadiationPattern.cxx)
target_link_libraries(dipoleRadiationPattern PUBLIC BasicFunctions FieldClasses ${ROOT_LIBRARIES} ${ROOTFFTW_LIB})

