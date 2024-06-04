cmake_minimum_required (VERSION 3.11)
project (mGEMS)
include(FetchContent)

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

## mGEMS library
add_library(libmgems
${CMAKE_CURRENT_SOURCE_DIR}/src/bin_reads.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/extract_bin.cpp
)
set_target_properties(libmgems PROPERTIES OUTPUT_NAME mgems)

## Dependencies
### OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  target_link_libraries(mGEMS OpenMP::OpenMP_CXX)
  target_link_libraries(libmgems OpenMP::OpenMP_CXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
endif()

### Check supported compression types
#### zlib
if ((DEFINED ZLIB_LIBRARY AND DEFINED ZLIB_INCLUDE_DIR) AND (NOT DEFINED ZLIB_FOUND))
  message(STATUS "zlib library provided in: " ${ZLIB_LIBRARY})
  message(STATUS "zlib headers provided in: " ${ZLIB_INCLUDE_DIR})
  include_directories(${ZLIB_INCLUDE_DIR})
  target_link_libraries(mGEMS ${ZLIB_LIBRARY})
  target_link_libraries(libmgems ${ZLIB_LIBRARY})
  set(MGEMS_HAVE_ZLIB 1)
else()
  find_package(ZLIB)
  if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIR})
    target_link_libraries(mGEMS ${ZLIB_LIBRARY})
    target_link_libraries(libmgems ${ZLIB_LIBRARY})
   set(MGEMS_HAVE_ZLIB 1)
  else()
    set(MGEMS_HAVE_ZLIB 0)
  endif()
endif()

#### bzip2
if (DEFINED BZIP2_LIBRARIES AND DEFINED BZIP2_INCLUDE_DIR AND (NOT DEFINED BZIP2_FOUND))
  message(STATUS "bzip2 library provided in: " ${BZIP2_LIBRARIES})
  message(STATUS "bzip2 headers provided in: " ${BZIP2_INCLUDE_DIR})
  include_directories(${BZIP2_INCLUDE_DIR})
  target_link_libraries(mGEMS ${BZIP2_LIBRARIES})
  target_link_libraries(libmgems ${BZIP2_LIBRARIES})
  set(MGEMS_HAVE_BZIP2 1)
else()
  find_package(BZip2)
  if (BZIP2_FOUND)
    include_directories(${BZIP2_INCLUDE_DIR})
    target_link_libraries(mGEMS ${BZIP2_LIBRARIES})
    target_link_libraries(libmgems ${BZIP2_LIBRARIES})
    set(MGEMS_HAVE_BZIP2 1)
  else()
    set(MGEMS_HAVE_BZIP2 0)
  endif()
endif()

#### lzma
if (DEFINED LIBLZMA_LIBRARY AND DEFINED LIBLZMA_INCLUDE_DIR AND (NOT DEFINED LIBLZMA_FOUND))
  message(STATUS "liblzma library provided in: " ${LIBLZMA_LIBRARY})
  message(STATUS "liblzma headers provided in: " ${LIBLZMA_INCLUDE_DIR})
  include_directories(${LIBLZMA_INCLUDE_DIR})
  target_link_libraries(mGEMS ${LIBLZMA_LIBRARY})
  target_link_libraries(libmgems ${LIBLZMA_LIBRARY})
  set(MGEMS_HAVE_LIBLZMA 1)
else()
  find_package(LibLZMA)
  if (LIBLZMA_FOUND)
    include_directories(${LIBLZMA_INCLUDE_DIR})
    target_link_libraries(mGEMS ${LIBLZMA_LIBRARY})
    target_link_libraries(libmgems ${LIBLZMA_LIBRARY})
    set(MGEMS_HAVE_LIBLZMA 1)
  else()
    set(MGEMS_HAVE_LIBLZMA 0)
  endif()
endif()

## bxzstr
if (DEFINED CMAKE_BXZSTR_HEADERS)
  message(STATUS "bxzstr headers provided in: ${CMAKE_BXZSTR_HEADERS}")
else()
  FetchContent_Declare(bxzstr
    GIT_REPOSITORY    https://github.com/tmaklin/bxzstr.git
    GIT_TAG           v1.1.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS        -D ZSTD_FOUND=0
                      -D ZLIB_FOUND=${MGEMS_HAVE_ZLIB}
		      -D BZIP2_FOUND=${MGEMS_HAVE_BZIP2}
		      -D LIBLZMA_FOUND=${MGEMS_HAVE_LIBLZMA}
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(bxzstr)
  set(CMAKE_BXZSTR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr/include)
endif()
include_directories(${CMAKE_BXZSTR_HEADERS})

## cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  message(STATUS "cxxargs headers provided in: ${CMAKE_CXXARGS_HEADERS}")
else()
  FetchContent_Declare(cxxargs
    GIT_REPOSITORY    https://github.com/tmaklin/cxxargs.git
    GIT_TAG           v1.1.4
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs"
    BUILD_IN_SOURCE   0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(cxxargs)
  set(CMAKE_CXXARGS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs/include)
endif()
include_directories("${CMAKE_CXXARGS_HEADERS}")

## Check dependencies and download them if not given
## cxxio
if (DEFINED CMAKE_CXXIO_HEADERS)
  message(STATUS "cxxio headers provided in: ${CMAKE_CXXIO_HEADERS}")
else()
  FetchContent_Declare(cxxio
    GIT_REPOSITORY    https://github.com/tmaklin/cxxio.git
    GIT_TAG           v0.1.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/cxxio"
    BUILD_IN_SOURCE   0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(cxxio)
  set(CMAKE_CXXIO_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxio/include)
endif()
include_directories("${CMAKE_CXXIO_HEADERS}")

## alignment-writer
if (DEFINED CMAKE_ALIGNMENT_WRITER_HEADERS AND DEFINED CMAKE_ALIGNMENT_WRITER_LIBRARY)
  message(STATUS "alignment-writer headers provided in: ${CMAKE_ALIGNMENT_WRITER_HEADERS}")
  message(STATUS "alignment-writer library provided in: ${CMAKE_ALIGNMENT_WRITER_LIBRARY}")
else()
  FetchContent_Declare(alignment-writer
    GIT_REPOSITORY    https://github.com/tmaklin/alignment-writer.git
    GIT_TAG           v0.5.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/alignment-writer"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_BXZSTR_HEADERS=${CMAKE_BXZSTR_HEADERS}
		    -D CMAKE_CXXARGS_HEADERS=${CMAKE_CXXARGS_HEADERS}
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(alignment-writer)
  add_dependencies(mGEMS libalignmentwriter)
  set_target_properties(alignment-writer PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/alignment-writer/include)
  set(CMAKE_ALIGNMENT_WRITER_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/lib/libalignment-writer.a)
endif()
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})
target_link_libraries(mGEMS ${CMAKE_ALIGNMENT_WRITER_LIBRARY})
target_link_libraries(libmgems ${CMAKE_ALIGNMENT_WRITER_LIBRARY})

## telescope
if (DEFINED CMAKE_TELESCOPE_LIBRARY AND DEFINED CMAKE_TELESCOPE_HEADERS)
  message(STATUS "telescope headers provided in: ${CMAKE_TELESCOPE_HEADERS}")
  message(STATUS "telescope library provided in: ${CMAKE_TELESCOPE_LIBRARY}")
else()
  FetchContent_Declare(telescope
    GIT_REPOSITORY    https://github.com/tmaklin/telescope.git
    GIT_TAG           v0.7.0-prerelease
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/telescope"
    BINARY_DIR        "${CMAKE_CURRENT_BINARY_DIR}/external/telescope"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS      -D CMAKE_BXZSTR_HEADERS=${CMAKE_BXZSTR_HEADERS}
		    -D CMAKE_CXXARGS_HEADERS=${CMAKE_CXXARGS_HEADERS}
		    -D CMAKE_CXXIO_HEADERS=${CMAKE_CXXIO_HEADERS}
		    -D CMAKE_ALIGNMENT_WRITER_HEADERS=${CMAKE_ALIGNMENT_WRITER_HEADERS}
		    -D CMAKE_ALIGNMENT_WRITER_LIBRARY=${CMAKE_ALIGNMENT_WRITER_LIBRARY}
		    -D CMAKE_BITMAGIC_HEADERS=${CMAKE_CURRENT_SOURCE_DIR}/external/telescope/external/BitMagic-7.12.3/src
		    -D CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		    -D "CMAKE_C_FLAGS=${CMAKE_C_FLAGS}"
		    -D "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
		    -D "CMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
		    -D "CMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(telescope)
  add_dependencies(mGEMS libalignmentwriter)
  add_dependencies(telescope libalignmentwriter)
  add_dependencies(mGEMS telescope)
  set_target_properties(telescope PROPERTIES EXCLUDE_FROM_ALL 1)
  set(CMAKE_TELESCOPE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/telescope/include)
  set(CMAKE_TELESCOPE_LIBRARY ${CMAKE_CURRENT_BINARY_DIR}/lib/libtelescope.a)
  set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/telescope/external/BitMagic-7.12.3/src)
endif()
include_directories(${CMAKE_TELESCOPE_HEADERS} ${CMAKE_BITMAGIC_HEADERS})
target_link_libraries(mGEMS ${CMAKE_TELESCOPE_LIBRARY})
target_link_libraries(libmgems ${CMAKE_TELESCOPE_LIBRARY})

## seamat
if (DEFINED CMAKE_SEAMAT_HEADERS)
  message(STATUS "seamat headers provided in: ${CMAKE_SEAMAT_HEADERS}")
else()
  FetchContent_Declare(seamat
    GIT_REPOSITORY    https://github.com/tmaklin/seamat.git
    GIT_TAG           v0.2.2
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/seamat"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS        -D CMAKE_BUILD_TESTS=0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(seamat)
  set(CMAKE_SEAMAT_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/seamat/include ${CMAKE_CURRENT_BINARY_DIR}/_deps/seamat-build/include)
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
  message( STATUS "mGEMS version: ${_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(_BUILD_VERSION "unknown")
endif()
string(TIMESTAMP _BUILD_TIMESTAMP)

## Generate a version.h file containing build version and timestamp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/version.h @ONLY)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link libraries
target_link_libraries(mGEMS libmgems)
