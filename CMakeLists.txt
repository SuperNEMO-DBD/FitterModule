# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(FitterModule)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
find_package(Falaise REQUIRED)

# Need ROOT library separately
find_package(ROOT REQUIRED)

# find the headers
include_directories(.)

# Build a dynamic library from our sources
add_library(Fitter SHARED fitter_module.h 
		   fitter_library.h 
		   fitter_library.cpp 
		   fitter_module.cpp)
target_link_libraries(Fitter PUBLIC Falaise::FalaiseModule)


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
add_executable(trial testing/trial.cpp)
target_link_libraries(trial PUBLIC Catch Fitter)

# CTest comand
add_test(NAME trial
  COMMAND trial -s
)
