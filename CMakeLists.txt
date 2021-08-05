project(opstop_cpp)
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
# set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS_DEBUG "-g3")
# project( gplines++ VERSION 0.0.1 DESCRIPTION "Generalized Splines in c++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(PYTHON_VERSION
    3.6
    CACHE STRING "Build bindings for Python version")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# include ifopt with custon ROS installation
set(ifopt_DIR /opt/ros/noetic/share/ifopt/cmake)
set(eigenpy_DIR /opt/ros/noetic/share/eigenpy/cmake)
set(pinocchio_DIR /opt/ros/noetic/share/cmake/pinocchio)

add_subdirectory(${CMAKE_SOURCE_DIR}/modules/gsplines_cpp)
find_package(GSL REQUIRED)
find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)

add_library(opstop SHARED ${CMAKE_SOURCE_DIR}/src/torque_constraint.cpp
                          ${CMAKE_SOURCE_DIR}/src/parametrization.cpp)

target_link_libraries(opstop ${GSL_LIBRARIES} ${Eigen3_LIBRARIES}
                      ${ifopt_LIBRARIES} ${pinocchio_LIBRARIES} gsplines)

target_include_directories(
  opstop
  PUBLIC ${CMAKE_SOURCE_DIR}/include
         ${GSL_INCLUDE_DIRS}
         ${Eigen3_INCLUDE_DIRS}
         ${ifopt_INCLUDE_DIRS}
         ${pinocchio_INCLUDE_DIRS}
         /usr/include/eigen3
         modules/gsplines_cpp/include/)

set_target_properties(opstop PROPERTIES LIBRARY_OUTPUT_DIRECTORY lib)

enable_testing()
add_subdirectory(tests)
