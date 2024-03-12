# project name and enable C++ support
project(rcbkdipole CXX C)

set(CMAKE_CXX_STANDARD 11)  # or 14, 17, etc.
# stop if cmake version is below 3.11
cmake_minimum_required(VERSION 3.11 FATAL_ERROR)

#find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development)

# Optimization on (for debug flags, change to Debug)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
include_directories("/opt/homebrew/anaconda3/include/python3.11")

add_subdirectory(src)
#add_subdirectory(pybind11)

# Find Pybind11 package
find_package(pybind11 REQUIRED)

# Add the Pybind11 include directories
include_directories(${pybind11_INCLUDE_DIRS})


pybind11_add_module(dipole_module src/dipole_amplitude.cpp)

# Set linker flags
target_link_libraries(dipole_module PRIVATE ${Python_LIBRARIES} amplitude)

set_target_properties(amplitude PROPERTIES
         LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
        ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib )

set_target_properties(dipole_module PROPERTIES
	 LIBRARY_OUTPUT_DIRECTORY ..
       	ARCHIVE_OUTPUT_DIRECTORY .. )

