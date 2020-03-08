cmake_minimum_required (VERSION 2.8.12)
project (mGEMS)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g -Wall -Wextra -Wpedantic")
endif()

## Set C++11 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 11)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++11)
endif()

## mGEMS executable
add_executable(mGEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

## Check supported compression types
find_package(BZip2)
if (BZIP2_FOUND)
  include_directories(${BZIP2_INCLUDE_DIRS})
  target_link_libraries(mGEMS ${BZIP2_LIBRARIES})
endif()    
find_package(LibLZMA)
if (LIBLZMA_FOUND)
  include_directories(${LIBLZMA_INCLUDE_DIRS})
  target_link_libraries(mGEMS ${LIBLZMA_LIBRARIES})
endif()
find_package(ZLIB)
if (ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(mGEMS ${ZLIB_LIBRARIES})
endif()

## Check dependencies and download them if not given
## telescope
if (DEFINED CMAKE_TELESCOPE_LIBRARY AND DEFINED CMAKE_TELESCOPE_HEADERS)
  find_library(TELESCOPE NAMES telescope HINTS ${CMAKE_TELESCOPE_LIBRARY})
  target_link_libraries(mGEMS ${TELESCOPE})
  include_directories("${CMAKE_TELESCOPE_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-telescope.txt.in ${CMAKE_BINARY_DIR}/external/telescope-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "CMake step for telescope failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "Build step for telescope failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/telescope
                   ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/build)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/telescope/include)
  set_target_properties(telescope PROPERTIES EXCLUDE_FROM_ALL 1)
  target_link_libraries(mGEMS libtelescope)
endif()

## bxzstr
if (DEFINED CMAKE_BXZSTR_HEADERS)
  include_directories("${CMAKE_BXZSTR_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-bxzstr.txt.in ${CMAKE_BINARY_DIR}/external/bxzstr-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "CMake step for bxzstr failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "Build step for bxzstr failed: ${result}")
  endif()
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr/include)
endif()

## cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  include_directories("${CMAKE_CXXARGS_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-cxxargs.txt.in ${CMAKE_BINARY_DIR}/external/cxxargs-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "CMake step for cxxargs failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "Build step for cxxargs failed: ${result}")
  endif()
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs/include)
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
configure_file(${CMAKE_SOURCE_DIR}/include/version.h.in ${CMAKE_BINARY_DIR}/include/version.h @ONLY)

## external/include has the version info
include_directories(include ${CMAKE_SOURCE_DIR}/external)

## mGEMS library
add_library(libmgems ${CMAKE_CURRENT_SOURCE_DIR}/src/mGEMS.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/bin_reads.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/extract_bin.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/mGEMS.cpp)
set_target_properties(libmgems PROPERTIES OUTPUT_NAME mgems)

# Link libraries
target_link_libraries(mGEMS libmgems)
