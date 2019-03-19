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
			      fitter_module.cpp)
target_link_libraries(Fitter PUBLIC Falaise::FalaiseModule)


# Build a dynamic library from our sources for algorithm only
# required for testing executables.
# include_directories(${ROOT_INCLUDE_DIRS})
# add_library(cdtorAlgorithm SHARED fitter_library.h 
#			  fitter_library.cpp)
# target_link_libraries(cdtorAlgorithm PUBLIC ${ROOT_LIBRARIES})


# Configure example pipeline script for use from the build dir
# configure_file("cdtorExample.conf.in" "fitterExample.conf" @ONLY)

# Build the testing code, tell CTest about it
# enable_testing()
# set(CMAKE_CXX_STANDARD 11)

# main catch library compile only once
# add_library(Catch testing/catch-main.cpp)
# target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/testing)

# add test executable
# include_directories(. testing)
# add_executable(trial testing/trial.cpp)
# target_link_libraries(trial PUBLIC Catch cdtorAlgorithm)

# CTest comand
# add_test(NAME trial
#   COMMAND trial -s
# )
