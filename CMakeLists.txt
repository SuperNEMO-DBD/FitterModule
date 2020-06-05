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
add_library(SNFitterModule SHARED fitter_module.cpp)
target_link_libraries(SNFitterModule PRIVATE SNFitter Falaise::FalaiseModule)

# Configure example pipeline script for use from the build dir
configure_file("fit.conf.in" "fit.conf" @ONLY)

# Build the testing code, tell CTest about it
enable_testing()

# main catch library compile only once
add_library(Catch testing/catch-main.cpp)
target_include_directories(Catch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/testing)

include_directories(${ROOT_INCLUDE_DIRS})

# add test executable
add_executable(linetest testing/linetest.cpp)
target_link_libraries(linetest PUBLIC Catch SNFitter ${ROOT_LIBRARIES})

# add test executable
add_executable(helixtest testing/helixtest.cpp)
target_link_libraries(helixtest PUBLIC Catch SNFitter ${ROOT_LIBRARIES})

# add test executable
add_executable(rhtest testing/rhtest.cpp)
target_link_libraries(rhtest PUBLIC Catch SNFitter ${ROOT_LIBRARIES})

# add test executable
add_executable(pftest testing/pftest.cpp)
target_link_libraries(pftest PUBLIC Catch SNFitter ${ROOT_LIBRARIES})

# add test executable
add_executable(blfittest testing/blfittest.cpp)
target_link_libraries(blfittest PUBLIC Catch SNFitter ${ROOT_LIBRARIES})

# add test executable
add_executable(linefiletest testing/linefiletest.cpp)
target_link_libraries(linefiletest PUBLIC Catch SNFitter ${ROOT_LIBRARIES})


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

