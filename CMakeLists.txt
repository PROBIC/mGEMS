cmake_minimum_required (VERSION 2.8.12)
project (msweep-assembly)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wno-sign-compare -g")

## Set C++11 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 11)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++11)
endif()

## Configure bxzstr
set(LIBLZMA_FOUND 0)
set(BZIP2_FOUND 0)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/bxzstr)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/telescope
                 ${CMAKE_BINARY_DIR}/external/telescope/build)

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY "${local_dir}"
    OUTPUT_VARIABLE _BUILD_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "mSWEEP-assembly plugin version: ${_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(_BUILD_VERSION "unknown")
endif()
string(TIMESTAMP _BUILD_TIMESTAMP)

## Generate a version.h file containing build version and timestamp
configure_file(${CMAKE_SOURCE_DIR}/include/version.h.in ${CMAKE_BINARY_DIR}/include/version.h @ONLY)

include_directories(include ${CMAKE_SOURCE_DIR}/external ${CMAKE_BINARY_DIR}/include)
file(GLOB ASSIGN_READS_SOURCES "${CMAKE_SOURCE_DIR}/src/assign_reads/*.cpp")
file(GLOB BUILD_SAMPLE_SOURCES "${CMAKE_SOURCE_DIR}/src/build_sample/*.cpp")
file(GLOB READ_ALIGNMENT_SOURCES "${CMAKE_SOURCE_DIR}/src/read_alignment/*.cpp")

add_executable(assign_reads ${ASSIGN_READS_SOURCES})
add_executable(build_sample ${BUILD_SAMPLE_SOURCES})
add_executable(read_alignment ${READ_ALIGNMENT_SOURCES})

# Link libraries
target_link_libraries(assign_reads z)
target_link_libraries(build_sample z)
target_link_libraries(read_alignment z)
