cmake_minimum_required(VERSION 3.18)
project(parRSB LANGUAGES C VERSION 0.9.3)

option(PARRSB_DEBUG "Build with debug symbols instead of -O2" OFF)
option(PARRSB_SHARED "Build a shared library instead of static" ON)
option(PARRSB_EXAMPLES "Build parRSB examples" ON)
option(PARRSB_TESTS "Build parRSB examples" ON)
set(PARRSB_SYNC "reduction" CACHE STRING "Sync mode: reduction | barrier")
set(PARRSB_VISIBILITY "hidden" CACHE STRING "Symbol visibility: hidden | default")

set(PARRSB_FLAGS
  -std=c99
  -Wall -Wextra -Wpedantic
  -Wno-unused-function
  -Wno-c23-extensions
  $<IF:$<BOOL:${PARRSB_DEBUG}>,-g -O0,-O2>
)

find_package(MPI REQUIRED)
find_package(gs REQUIRED)

#Link:https://gitlab.kitware.com/cmake/community/wikis/doc/cmake/RPATH-handling
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/lib;${gs_LIBRARY_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


file(GLOB SRCS  "${CMAKE_SOURCE_DIR}/src/*.c")
if(PARRSB_SHARED)
  add_library(parRSB SHARED ${SRCS})
else()
  add_library(parRSB STATIC ${SRCS})
endif()
set_target_properties(parRSB PROPERTIES OUTPUT_NAME parRSB)
set_target_properties(parRSB PROPERTIES VERSION ${PROJECT_VERSION})
target_include_directories(parRSB PUBLIC $<INSTALL_INTERFACE:include>
  PRIVATE $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_options(parRSB PRIVATE ${PARRSB_FLAGS})
target_link_libraries(parRSB PRIVATE gs::gs MPI::MPI_C m)

if(PARRSB_SYNC STREQUAL "reduction")
  target_compile_definitions(parRSB PRIVATE PARRSB_SYNC_REDUCTION)
elseif(PARRSB_SYNC STREQUAL "barrier")
  target_compile_definitions(parRSB PRIVATE PARRSB_SYNC_BARRIER)
endif()

if(NOT PARRSB_VISIBILITY STREQUAL "")
  target_compile_definitions(parRSB PRIVATE PARRSB_VISIBILITY=${PARRSB_VISIBILITY})
endif()

function(add_binary src)
  get_filename_component(name "${src}" NAME_WE)
  add_executable(${name} "${src}")
  target_compile_options(${name} PRIVATE ${PARRSB_FLAGS})
  target_link_libraries(${name} PRIVATE parRSB gs::gs MPI::MPI_C m)
  install(TARGETS ${name} RUNTIME DESTINATION bin)
  target_include_directories(${name} PRIVATE "${CMAKE_SOURCE_DIR}/src")
endfunction()

if(PARRSB_EXAMPLES)
  file(GLOB EXAMPLE_SRCS "${CMAKE_SOURCE_DIR}/example/*.c")
  foreach(src ${EXAMPLE_SRCS})
    add_binary(${src})
  endforeach()
endif()

if(PARRSB_EXAMPLES)
  file(GLOB EXAMPLE_SRCS "${CMAKE_SOURCE_DIR}/test/*.c")
  foreach(src ${EXAMPLE_SRCS})
    add_binary(${src})
  endforeach()
endif()

# install parRSB
install(TARGETS parRSB
  EXPORT parRSBTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
file(GLOB HEADERS "${CMAKE_SOURCE_DIR}/src/*.h")
install(FILES ${HEADERS} DESTINATION include)

# export targets/configs
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

# `format` target
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
