##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================
cmake_minimum_required(VERSION 3.8...3.15 FATAL_ERROR)
project(advection CXX)


list(APPEND CMAKE_MODULE_PATH "~/repositories/CMakeFindModules")

# Find VAPOR packages
# https://github.com/ayenpure/CMakeFindModules/blob/master/FindVAPOR.cmake
find_package(VAPOR REQUIRED QUIET)
# https://github.com/ayenpure/CMakeFindModules/blob/master/FindVAPORDEPS.cmake
find_package(VAPORDEPS REQUIRED QUIET)
# Find Boost
find_package(Boost COMPONENTS program_options filesystem REQUIRED)
# Find VTK-m for comparison of output field
find_package(VTK REQUIRED)
find_package(VTKm REQUIRED QUIET)

# OMP optional on others, requred locally
find_package(OMP REQUIRED QUIET)

add_executable(benchmark benchmark.cxx Advection.h Advection.cpp)
set_target_properties(benchmark PROPERTIES CXX_STANDARD 11)
target_include_directories(benchmark PRIVATE ${VAPORDEPS_INCLUDE_DIRS} ${VAPOR_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${OMP_INCLUDE})
target_link_libraries(benchmark PRIVATE ${VAPOR_LIBRARIES} ${Boost_LIBRARIES} ${OMP_LIBRARY})

add_executable(ftle ftle.cxx FTLEHelper.h Advection.h Advection.cpp)
set_target_properties(ftle PROPERTIES CXX_STANDARD 11)
target_include_directories(ftle PRIVATE ${VAPORDEPS_INCLUDE_DIRS} ${VAPOR_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
target_link_libraries(ftle PRIVATE ${VAPOR_LIBRARIES} ${Boost_LIBRARIES})

include(${VTK_USE_FILE})
add_executable(compare CompareWithVTKm.cxx Advection.h Advection.cpp)
set_target_properties(compare PROPERTIES CXX_STANDARD 11)
target_include_directories(compare PRIVATE ${VAPORDEPS_INCLUDE_DIRS} ${VAPOR_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
target_link_libraries(compare PRIVATE ${VAPOR_LIBRARIES} ${Boost_LIBRARIES} ${VTK_LIBRARIES} vtkm_cont vtkm_io vtkm_worklet)

#add_executable(convert ConvertToVTK.cxx)
#set_target_properties(convert PROPERTIES CXX_STANDARD 11)
#target_include_directories(convert PRIVATE ${VAPORDEPS_INCLUDE_DIRS} ${VAPOR_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
#target_link_libraries(convert PRIVATE ${VAPOR_LIBRARIES} ${Boost_LIBRARIES})
