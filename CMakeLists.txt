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

# Find the VTK-m package
find_package(VAPOR REQUIRED QUIET)
find_package(OpenGL REQUIRED)
INCLUDE_DIRECTORIES(${VAPOR_INCLUDE_DIRS} ${OpenGL_INCLUDE_DIRS})

add_executable(advection advection.cxx)
target_link_libraries(advection PRIVATE ${VAPOR_LIBRARIES} ${OpenGL_LIBRARIES})

