cmake_minimum_required (VERSION 2.8.12)
project (mGEMS)

set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
elseif(CMAKE_BUILD_TYPE MATCHES Debug)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -Wall -Wextra -Wpedantic")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g -Wall -Wextra -Wpedantic")
endif()

## Set C++17 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 17)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++17)
endif()

## mGEMS executable
add_executable(mGEMS ${CMAKE_CURRENT_SOURCE_DIR}/src/mGEMS.cpp)

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

## bxzstr
if (DEFINED CMAKE_BXZSTR_HEADERS)
  message(STATUS "bxzstr headers provided in: ${CMAKE_BXZSTR_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-bxzstr.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "CMake step for bxzstr failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "Build step for bxzstr failed: ${result}")
  endif()
  set(CMAKE_BXZSTR_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr/include)
endif()
include_directories(${CMAKE_BXZSTR_HEADERS})

## cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  message(STATUS "cxxargs headers provided in: ${CMAKE_CXXARGS_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-cxxargs.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "CMake step for cxxargs failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "Build step for cxxargs failed: ${result}")
  endif()
  set(CMAKE_CXXARGS_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs/include)
endif()
include_directories("${CMAKE_CXXARGS_HEADERS}")

## cxxio
if (DEFINED CMAKE_CXXIO_HEADERS)
  message(STATUS "cxxio headers provided in: ${CMAKE_CXXIO_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-cxxio.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio-download )
  if(result)
    message(FATAL_ERROR "CMake step for cxxio failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio-download )
  if(result)
    message(FATAL_ERROR "Build step for cxxio failed: ${result}")
  endif()
  set(CMAKE_CXXIO_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/cxxio/include)
endif()
include_directories("${CMAKE_CXXIO_HEADERS}")

## alignment-writer
if (DEFINED CMAKE_ALIGNMENT_WRITER_HEADERS)
  message(STATUS "alignment-writer headers provided in: ${CMAKE_ALIGNMENT_WRITER_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-alignment-writer.txt.in ${CMAKE_BINARY_DIR}/external/alignment-writer-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/alignment-writer-download )
  if(result)
    message(FATAL_ERROR "CMake step for alignment-writer failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/alignment-writer-download )
  if(result)
    message(FATAL_ERROR "Build step for alignment-writer failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer
                   ${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer/build)
  set_target_properties(alignment-writer PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer/include)
endif()
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})

## Check dependencies and download them if not given
## telescope
if (DEFINED CMAKE_TELESCOPE_HEADERS)
  target_link_libraries(mGEMS libtelescope)
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-telescope.txt.in ${CMAKE_CURRENT_BINARY_DIR}/external/telescope-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "CMake step for telescope failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/telescope-download )
  if(result)
    message(FATAL_ERROR "Build step for telescope failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/telescope
                   ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/build)
  set_target_properties(telescope PROPERTIES EXCLUDE_FROM_ALL 1)
  target_link_libraries(mGEMS libtelescope libalignmentwriter)
  set(CMAKE_TELESCOPE_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/include ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/build/include)
  if (NOT DEFINED CMAKE_BITMAGIC_HEADERS)
    set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/telescope/external/BitMagic-7.12.3/src)
  endif()
endif()

include_directories(${CMAKE_TELESCOPE_HEADERS} ${CMAKE_BITMAGIC_HEADERS})

## seamat
if (DEFINED CMAKE_SEAMAT_HEADERS)
  message(STATUS "seamat headers provided in: ${CMAKE_SEAMAT_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-seamat.txt.in ${CMAKE_BINARY_DIR}/external/seamat-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/seamat-download )
  if(result)
    message(FATAL_ERROR "CMake step for seamat failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/seamat-download )
  if(result)
    message(FATAL_ERROR "Build step for seamat failed: ${result}")
  endif()
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/external/seamat
                   ${CMAKE_CURRENT_BINARY_DIR}/external/seamat/build)
  target_link_libraries(mGEMS)
  set(CMAKE_SEAMAT_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/seamat/include)
endif()
include_directories(${CMAKE_SEAMAT_HEADERS})

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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/version.h @ONLY)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

## mGEMS library
add_library(libmgems
${CMAKE_CURRENT_SOURCE_DIR}/src/bin_reads.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/extract_bin.cpp
)
set_target_properties(libmgems PROPERTIES OUTPUT_NAME mgems)

# Link libraries
target_link_libraries(mGEMS libmgems)
