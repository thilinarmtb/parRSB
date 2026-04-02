cmake_minimum_required(VERSION 3.18)
project(parRSB LANGUAGES C VERSION 0.9.3)

option(PARRSB_DEBUG "Build with debug symbols instead of -O2" OFF)
option(PARRSB_SHARED "Build a shared library instead of static" ON)
set(PARRSB_SYNC "reduction" CACHE STRING "Sync mode: reduction | barrier")
set(PARRSB_VISIBILITY "hidden" CACHE STRING "Symbol visibility: hidden | default")

find_package(MPI REQUIRED)
find_package(gs REQUIRED)

file(GLOB SRCS  "${CMAKE_SOURCE_DIR}/src/*.c")
if(PARRSB_SHARED)
  add_library(parRSB SHARED ${SRCS})
else()
  add_library(parRSB STATIC ${SRCS})
endif()
set_target_properties(parRSB PROPERTIES OUTPUT_NAME parRSB)
set_target_properties(parRSB PROPERTIES VERSION ${PROJECT_VERSION}
    PUBLIC_HEADER src/parrsb.h)
target_include_directories(parRSB PRIVATE "${CMAKE_SOURCE_DIR}/src")
target_compile_options(parRSB PRIVATE
  -std=c99
  -Wall -Wextra -Wpedantic
  -Wno-unused-function
  -Wno-c23-extensions
  $<IF:$<BOOL:${PARRSB_DEBUG}>,-g -O0,-O2>
)
target_link_libraries(parRSB PRIVATE gs::gs MPI::MPI_C m)

if(PARRSB_SYNC STREQUAL "reduction")
  target_compile_definitions(parRSB PRIVATE PARRSB_SYNC_REDUCTION)
elseif(PARRSB_SYNC STREQUAL "barrier")
  target_compile_definitions(parRSB PRIVATE PARRSB_SYNC_BARRIER)
endif()

if(NOT PARRSB_VISIBILITY STREQUAL "")
  target_compile_definitions(parRSB PRIVATE PARRSB_VISIBILITY=${PARRSB_VISIBILITY})
endif()

function(add_bin source_file)
  get_filename_component(name "${source_file}" NAME_WE)
  add_executable(${name} "${source_file}")
  target_compile_options(${name} PRIVATE
    -std=c99
    -Wall -Wextra -Wpedantic
    -Wno-unused-function
    -Wno-c23-extensions
    $<IF:$<BOOL:${PARRSB_DEBUG}>,-g -O0,-O2>
  )
  if(PARRSB_SHARED)
    target_link_libraries(${name} PRIVATE parRSB gs::gs MPI::MPI_C m)
  else()
    target_link_libraries(${name} PRIVATE parRSB gs::gs MPI::MPI_C m)
  endif()
  target_include_directories(${name} PRIVATE "${CMAKE_SOURCE_DIR}/src")
  add_dependencies(${name} parRSB)
endfunction()

file(GLOB EXAMPLE_SRCS "${CMAKE_SOURCE_DIR}/example/*.c")
foreach(src ${EXAMPLE_SRCS})
  add_bin(${src})
endforeach()

file(GLOB TEST_SRCS "${CMAKE_SOURCE_DIR}/test/*.c")
foreach(src ${TEST_SRCS})
  add_bin(${src})
endforeach()

install(TARGETS parRSB
  EXPORT parRSBTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
)

set(EXPORT_TARGETS_FILE "parRSBTargets.cmake")
install(EXPORT parRSBTargets
  FILE parRSBTargets.cmake
  NAMESPACE parRSB::
  DESTINATION lib/cmake/parRSB
)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/parRSBConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
configure_file(
  "${CMAKE_SOURCE_DIR}/parRSBConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/parRSBConfig.cmake"
  @ONLY
)
install(FILES
  "${CMAKE_BINARY_DIR}/parRSBConfig.cmake"
  "${CMAKE_BINARY_DIR}/parRSBConfigVersion.cmake"
  DESTINATION lib/cmake/parRSB
)

find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
  file(GLOB_RECURSE ALL_C_FILES
    "${CMAKE_SOURCE_DIR}/src/*.c"
    "${CMAKE_SOURCE_DIR}/src/*.h"
    "${CMAKE_SOURCE_DIR}/example/*.c"
    "${CMAKE_SOURCE_DIR}/test/*.c"
  )
  add_custom_target(format
    COMMAND ${CLANG_FORMAT} -i ${ALL_C_FILES}
    COMMENT "Running clang-format"
  )
endif()
