include(GoogleTest)
link_directories(${PROJECT_SOURCE_DIR}/build)
link_libraries(${GSL_LIBRARIES} ${Eigen3_LIBRARIES} gsplines opstop gtest
               pthread)
include_directories(PUBLIC ${PROJECT_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS}
                    ${Eigen3_LIBRARIES})

file(GLOB_RECURSE test_list ${PROJECT_SOURCE_DIR}/tests/*.cpp)
file(GLOB_RECURSE test_py_list ${PROJECT_SOURCE_DIR}/tests/*.py)
list(REMOVE_ITEM test_py_list "${PROJECT_SOURCE_DIR}/tests/__init__.py")
list(REMOVE_ITEM test_py_list "${PROJECT_SOURCE_DIR}/tests/tools.py")

foreach(file_path ${test_list})
  get_filename_component(barename ${file_path} NAME)
  string(REPLACE ".cpp" "" new_name ${barename})
  add_executable(${new_name} ${file_path})
  add_test(
    NAME ${new_name}
    COMMAND ${new_name}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
endforeach()

foreach(file_path ${test_py_list})
  get_filename_component(barename ${file_path} NAME)
  string(REPLACE ".py" "_py" new_name ${barename})
  string(REPLACE ".py" "" mod_name ${barename})
  add_test(
    NAME ${new_name}
    COMMAND python3 -m unittest tests.${mod_name}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})

endforeach()

add_executable(gtest_l2_constraints
               ${PROJECT_SOURCE_DIR}/tests/jerk_l2_constraints.cpp)
gtest_discover_tests(gtest_l2_constraints)
add_executable(gtest_minimum_time_bounded_jerk_l2
               ${PROJECT_SOURCE_DIR}/tests/minimum_time_bounded_jerk_l2.cpp)
gtest_discover_tests(gtest_minimum_time_bounded_jerk_l2
                     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/tests)
