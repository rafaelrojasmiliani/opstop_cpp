project(opstop_cpp)
include(CTest)
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

add_subdirectory(${PROJECT_SOURCE_DIR}/modules/pybind11)
add_subdirectory(${PROJECT_SOURCE_DIR}/modules/gsplines_cpp)
find_package(GSL REQUIRED)
find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)

add_library(
  opstop SHARED
  ${PROJECT_SOURCE_DIR}/src/torque_constraint.cpp
  ${PROJECT_SOURCE_DIR}/src/parametrization.cpp
  ${PROJECT_SOURCE_DIR}/src/time_cost.cpp
  ${PROJECT_SOURCE_DIR}/src/parametrization_variables.cpp
  ${PROJECT_SOURCE_DIR}/src/velocity_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/acceleration_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/jerk_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/diffeo_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/ipopt_problem.cpp
  ${PROJECT_SOURCE_DIR}/src/differ.cpp)

target_link_libraries(opstop ${GSL_LIBRARIES} ${Eigen3_LIBRARIES}
                      ${ifopt_LIBRARIES} ${pinocchio_LIBRARIES} gsplines)

target_include_directories(
  opstop
  PUBLIC ${PROJECT_SOURCE_DIR}/include
         ${GSL_INCLUDE_DIRS}
         ${Eigen3_INCLUDE_DIRS}
         ${ifopt_INCLUDE_DIRS}
         ${pinocchio_INCLUDE_DIRS}
         /usr/include/eigen3
         modules/gsplines_cpp/include/)

set_target_properties(opstop PROPERTIES LIBRARY_OUTPUT_DIRECTORY lib)

pybind11_add_module(pyopstop ${PROJECT_SOURCE_DIR}/bindings/bindings.cpp)
target_include_directories(pyopstop PUBLIC include /usr/include/eigen3
                                           modules/gsplines_cpp/include/)
target_link_libraries(pyopstop PUBLIC ${Eigen3_LIBRARIES} ${ifopt_LIBRARIES}
                                      gsplines opstop)

set_target_properties(pyopstop PROPERTIES OUTPUT_NAME opstop)
if(${PROJECT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  enable_testing()
  add_subdirectory(tests)
endif()

if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
  execute_process(
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
      ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)
endif()
