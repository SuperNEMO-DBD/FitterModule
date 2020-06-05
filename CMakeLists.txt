# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 3.3)
project(FitterModule)

# Build settings
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Dependencies
find_package(Falaise REQUIRED)

# Falaise will expose ROOT for us, but ensure we
# have the libraries we need.
find_package(ROOT REQUIRED Spectrum Minuit2 Minuit)

# Build a dynamic library from our sources
add_library(SNFitter SHARED
  SNFitter/DataProducts.h
  SNFitter/SNFitter.cpp
  SNFitter/SNFitter.h
  SNFitter/SNFitterImpl.cpp
  SNFitter/SNFitterImpl.h)
target_include_directories(SNFitter
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(SNFitter PRIVATE ${ROOT_LIBRARIES})

# Build the Falaise plugin
add_library(SNFitterModule SHARED SNFitterModule.cpp)
target_link_libraries(SNFitterModule PRIVATE SNFitter Falaise::FalaiseModule)

# Configure example pipeline script for use from the build dir
configure_file("fit.conf.in" "fit.conf" @ONLY)

# Build the testing code, tell CTest about it
enable_testing()

# main catch library compile only once
add_library(Catch testing/catch-main.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/testing)

include_directories(${ROOT_INCLUDE_DIRS})

foreach(__testname linetest linefiletest helixtest rhtest pftest blfittest)
  add_executable(${__testname} testing/${__testname}.cpp)
  target_link_libraries(${__testname} PRIVATE Catch SNFitter ${ROOT_LIBRARIES})

  add_test(NAME ${__testname} COMMAND ${__testname} -s)
endforeach()
