cmake_minimum_required(VERSION 3.5)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

project("SLHD")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(SLHD_lib
  "src/design.cpp"
  "src/design.h"
  "src/utils.h"
  "src/SLHD.cpp"
  "src/SLHD.h"
)

add_library(manager_lib
  "src/design.cpp"
  "src/design.h"
  "src/utils.h"
  "src/SLHD.cpp"
  "src/SLHD.h"
  "src/manager.cpp"
  "src/manager.h"
)

include_directories(${PROJECT_SOURCE_DIR}/src)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/../bin)

add_executable(SLHD_test
  "src/SLHD_test.cpp"
)

add_executable(manager_test
  "src/manager_test.cpp"
)

target_link_libraries(SLHD_test SLHD_lib)

target_link_libraries(manager_test manager_lib)