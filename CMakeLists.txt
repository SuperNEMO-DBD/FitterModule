# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(FitterModule)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Falaise REQUIRED)

# Need ROOT library separately
find_package(ROOT REQUIRED Spectrum Minuit2 Minuit)

# find the headers
include_directories(.)

# Build a dynamic library from our sources
add_library(Fitter SHARED fitter_module.h 
		   fitter_library.h 
		   fitter_library.cpp 
		   fitter_module.cpp)
target_link_libraries(Fitter PUBLIC Falaise::FalaiseModule ${ROOT_LIBRARIES})

# Build a dynamic library from our sources for algorithm only
# required for testing executables.
include_directories(${ROOT_INCLUDE_DIRS})
add_library(FitAlgorithm SHARED fitter_library.h 
			  fitter_library.cpp)
target_link_libraries(FitAlgorithm PUBLIC ${ROOT_LIBRARIES})


# Configure example pipeline script for use from the build dir
configure_file("fit.conf.in" "fit.conf" @ONLY)

# Build the testing code, tell CTest about it
enable_testing()
set(CMAKE_CXX_STANDARD 11)

# main catch library compile only once
add_library(Catch testing/catch-main.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/testing)

# add test executable
include_directories(. testing)
add_executable(linetest testing/linetest.cpp)
target_link_libraries(linetest PUBLIC Catch FitAlgorithm)

# add test executable
add_executable(helixtest testing/helixtest.cpp)
target_link_libraries(helixtest PUBLIC Catch FitAlgorithm)

# add test executable
add_executable(rhtest testing/rhtest.cpp)
target_link_libraries(rhtest PUBLIC Catch FitAlgorithm)

# add test executable
add_executable(pftest testing/pftest.cpp)
target_link_libraries(pftest PUBLIC Catch FitAlgorithm)

# add test executable
add_executable(blfittest testing/blfittest.cpp)
target_link_libraries(blfittest PUBLIC Catch FitAlgorithm)

# add test executable
add_executable(linefiletest testing/linefiletest.cpp)
target_link_libraries(linefiletest PUBLIC Catch FitAlgorithm)


# CTest comand
add_test(NAME linetest
  COMMAND linetest -s
)

# CTest comand
add_test(NAME rhtest
  COMMAND rhtest -s
)

# CTest comand
add_test(NAME pftest
  COMMAND pftest -s
)

# CTest comand
add_test(NAME helixtest
  COMMAND helixtest -s
)

# CTest comand
add_test(NAME linefiletest
  COMMAND linefiletest -s
)

# CTest comand
add_test(NAME blfittest
  COMMAND blfittest -s
)

