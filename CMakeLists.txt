cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(opstop VERSION "0.0.1")
include(CTest)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops -mfpmath=sse -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -pthread ")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_INSTALL_PREFIX /usr)
set(PYBIND11_FINDPYTHON
    OFF
    CACHE BOOL "")

set(PYBIND11_PYTHON_VERSION
    "3.8"
    CACHE STRING "")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT TARGET pybind11)
  add_subdirectory(${PROJECT_SOURCE_DIR}/modules/pybind11)
endif()
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

target_link_libraries(opstop ${EIGEN3_LIBRARIES} ${ifopt_LIBRARIES}
                      ${pinocchio_LIBRARIES} ${gsplines_LIBRARIES})

target_include_directories(
  opstop
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>
         ${EIGEN3_INCLUDE_DIRS}
         ${ifopt_INCLUDE_DIRS}
         ${pinocchio_INCLUDE_DIRS}
         ${gsplines_INCLUDE_DIRS})

set_target_properties(
  opstop
  PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON
             VERSION ${opstop_VERSION}
             SOVERSION ${opstop_VERSION_MAJOR}.${opstop_VERSION_MINOR}
             EXPORT_NAME OpStop)

pybind11_add_module(pyopstop ${PROJECT_SOURCE_DIR}/bindings/bindings.cpp)

target_include_directories(
  pyopstop
  PUBLIC ${EIGEN3_INCLUDE_DIRS}
         ${gsplines_INCLUDE_DIRS}
         $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/modules/pybind11>
         $<INSTALL_INTERFACE:include>)

target_link_libraries(
  pyopstop PUBLIC ${EIGEN3_LIBRARIES} ${ifopt_LIBRARIES} ${pinocchio_LIBRARIES}
                  ${gsplines_LIBRARIES} opstop)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR}
   AND ${BUILD_TESTING})
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

include(GNUInstallDirs)
set(INSTALL_CMAKE_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/OpStop)

install(
  TARGETS opstop
  EXPORT OpStopTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import site; print(site.getsitepackages()[0])"
  OUTPUT_VARIABLE PYTHON_INSTALL_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

install(DIRECTORY python/ DESTINATION ${PYTHON_INSTALL_PATH})
install(TARGETS pyopstop DESTINATION ${PYTHON_INSTALL_PATH}/opstop)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(EXPORT OpStopTargets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/OpStopTargets.cmake)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/OpStopConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)

configure_package_config_file(
  cmake/OpStopConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/OpStopConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CMAKE_CONFIG_DIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/OpStopConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/OpStopConfigVersion.cmake
        DESTINATION ${INSTALL_CMAKE_CONFIG_DIR})

install(EXPORT OpStopTargets DESTINATION ${INSTALL_CMAKE_CONFIG_DIR})

# Packaging
set(CPACK_PACKAGE_VENDOR "Rafael A. Rojas")
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VERSION ${opstop_VERSION})
set(CPACK_SYSTEM_NAME ${CMAKE_HOST_SYSTEM_PROCESSOR})

# Debian
find_program(DPKG_PROG dpkg DOC "'dpkg' executable")
if(DPKG_PROG)
  execute_process(
    COMMAND ${DPKG_PROG} --print-architecture
    OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  # Change system name to use the correct architecture in file name
  set(CPACK_SYSTEM_NAME ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE})
endif()
# Debian versions require a dash
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION}-1)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Rafael A. Rojas")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "coinor-libipopt-dev, ros-noetic-ifopt, libeigen3-dev, ros-noetic-hpp-fcl, robotpkg-pinocchio, gsplines"
)

include(CPack)
# set(PYTHON_EXAMPLE_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/examuexample.py)
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/example_template.py
# ${PYTHON_EXAMPLE_SCRIPT} @ONLY) file(WRITE ${PYTHON_EXAMPLE_SCRIPT} "#!
# /usr/bin/env python\nimport
# sys\nsys.path.append('${CMAKE_LIBRARY_OUTPUT_DIRECTORY}')\n") file(APPEND
# ${PYTHON_EXAMPLE_SCRIPT} "import your_pybind_module\n")
# execute_process(COMMAND chmod +x ${PYTHON_EXAMPLE_SCRIPT})
