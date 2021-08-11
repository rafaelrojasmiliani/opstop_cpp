link_directories(${CMAKE_SOURCE_DIR}/build)
link_libraries(${GSL_LIBRARIES} ${Eigen3_LIBRARIES} gsplines opstop)
include_directories(PUBLIC ${CMAKE_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS}
                    ${Eigen3_LIBRARIES})

file(GLOB_RECURSE test_list ${CMAKE_SOURCE_DIR}/tests/*.cpp)
file(GLOB_RECURSE test_py_list ${CMAKE_SOURCE_DIR}/tests/*.py)
list(REMOVE_ITEM test_py_list "__init__.py")

foreach(file_path ${test_list})
  get_filename_component(barename ${file_path} NAME)
  string(REPLACE ".cpp" "" new_name ${barename})
  add_executable(${new_name} ${file_path})
  add_test(
    NAME ${new_name}
    COMMAND ${new_name}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
endforeach()

foreach(file_path ${test_py_list})
  get_filename_component(barename ${file_path} NAME)
  string(REPLACE ".py" "" new_name ${barename})
  add_test(
    NAME ${new_name}
    COMMAND python3 -m unittest ${NEw_name}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests)
endforeach()
