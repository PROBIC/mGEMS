cmake_minimum_required (VERSION 2.8.12)
project (assign_reads)

## Set C++11 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 11)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++11)
endif()

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
configure_file(include/version.h.in ${CMAKE_SOURCE_DIR}/include/version.h @ONLY)

include_directories(include)
file(GLOB ASSIGN_READS_SOURCES "src/assign_reads/*.cpp")
file(GLOB BUILD_SAMPLE_SOURCES "src/build_sample/*.cpp")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

add_executable(assign_reads ${ASSIGN_READS_SOURCES})
add_executable(build_sample ${BUILD_SAMPLE_SOURCES})

# Link libraries
target_link_libraries(assign_reads z)
target_link_libraries(build_sample z)

# Configure compiler flags, release enables -O2 optimization.
set(CMAKE_BUILD_TYPE Release)
