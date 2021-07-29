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

add_subdirectory(${CMAKE_SOURCE_DIR}/modules/gsplines_cpp)
find_package(GSL REQUIRED)
find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(opstop SHARED ${CMAKE_SOURCE_DIR}/src/parametrization.cpp)

target_link_libraries(opstop ${GSL_LIBRARIES} ${Eigen3_LIBRARIES}
                      ${ifopt_LIBRARIES} gsplines)

target_include_directories(
  opstop
  PUBLIC ${CMAKE_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS}
         ${ifopt_INCLUDE_DIRS} /usr/include/eigen3
         modules/gsplines_cpp/include/)

set_target_properties(opstop PROPERTIES LIBRARY_OUTPUT_DIRECTORY lib)

# set(PYBIND11_PYTHON_VERSION "3.6") link_directories(${CMAKE_SOURCE_DIR}/build)
# add_subdirectory(${CMAKE_SOURCE_DIR}/modules/pybind11)
# pybind11_add_module(pyopstop ${CMAKE_SOURCE_DIR}/bindings/bindings.cpp)
# target_include_directories(pyopstop PUBLIC include ${GSL_INCLUDE_DIRS}
# /usr/include/eigen3) target_link_libraries(pyopstop PUBLIC ${GSL_LIBRARIES}
# ${Eigen3_LIBRARIES} ${ifopt_LIBRARIES} gsplines) enable_testing()
# add_subdirectory(tests)
