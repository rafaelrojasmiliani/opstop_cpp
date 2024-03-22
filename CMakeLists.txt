project(opstop DESCRIPTION "Time-optimal stop along gsplines")
include(CTest)
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE
    "-O3 -DNDEBUG -funroll-loops   -mfpmath=sse  -fopenmp")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(PYTHON_VERSION
    3.8
    CACHE STRING "Build bindings for Python version")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT TARGET pybind11)
  set(PYBIND11_PYTHON_VERSION "3.8")
  add_subdirectory(${PROJECT_SOURCE_DIR}/modules/pybind11)
endif()
find_package(GSL REQUIRED)
find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)
find_package(GSplines REQUIRED)

add_library(
  opstop SHARED
  ${PROJECT_SOURCE_DIR}/src/torque_constraint.cpp
  ${PROJECT_SOURCE_DIR}/src/parametrization.cpp
  ${PROJECT_SOURCE_DIR}/src/time_cost.cpp
  ${PROJECT_SOURCE_DIR}/src/excursion_cost.cpp
  ${PROJECT_SOURCE_DIR}/src/parametrization_variables.cpp
  ${PROJECT_SOURCE_DIR}/src/velocity_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/acceleration_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/jerk_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/jerk_l2_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/diffeo_constraints.cpp
  ${PROJECT_SOURCE_DIR}/src/ipopt_problem.cpp
  ${PROJECT_SOURCE_DIR}/src/differ.cpp)

target_link_libraries(
  opstop ${GSL_LIBRARIES} ${Eigen3_LIBRARIES} ${ifopt_LIBRARIES}
  ${pinocchio_LIBRARIES} ${gsplines_LIBRARIES})

target_include_directories(
  opstop
  PUBLIC ${PROJECT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS}
         ${EIGEN3_INCLUDE_DIRS} ${ifopt_INCLUDE_DIRS} ${pinocchio_INCLUDE_DIRS}
         ${gsplines_INCLUDE_DIRS})

set_target_properties(opstop PROPERTIES LIBRARY_OUTPUT_DIRECTORY lib)

pybind11_add_module(pyopstop ${PROJECT_SOURCE_DIR}/bindings/bindings.cpp)
target_include_directories(pyopstop PUBLIC include /usr/include/eigen3
                                           ${gsplines_INCLUDE_DIRS})

target_link_libraries(
  pyopstop PUBLIC ${Eigen3_LIBRARIES} ${ifopt_LIBRARIES} ${pinocchio_LIBRARIES}
                  ${gsplines_LIBRARIES} opstop)

set_target_properties(pyopstop PROPERTIES OUTPUT_NAME opstop)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
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

# set(PYTHON_EXAMPLE_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/examuexample.py)
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/example_template.py
# ${PYTHON_EXAMPLE_SCRIPT} @ONLY) file(WRITE ${PYTHON_EXAMPLE_SCRIPT} "#!
# /usr/bin/env python\nimport
# sys\nsys.path.append('${CMAKE_LIBRARY_OUTPUT_DIRECTORY}')\n") file(APPEND
# ${PYTHON_EXAMPLE_SCRIPT} "import your_pybind_module\n")
# execute_process(COMMAND chmod +x ${PYTHON_EXAMPLE_SCRIPT})
