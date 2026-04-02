cmake_minimum_required(VERSION 3.18)
project(parRSB C)

# Options
option(DEBUG "Build with debug symbols instead of -O2" OFF)
set(PARRSB_SYNC "reduction" CACHE STRING "Sync mode: reduction | barrier | (empty)")
set_property(CACHE PARRSB_SYNC PROPERTY STRINGS reduction barrier "")
option(PARRSB_SHARED "Build a shared library instead of static" OFF)
set(PARRSB_VISIBILITY "hidden" CACHE STRING "Symbol visibility: hidden | default | (empty)")
set_property(CACHE PARRSB_VISIBILITY PROPERTY STRINGS hidden default "")

# MPI (replaces the bare `mpicc` compiler driver)
find_package(MPI REQUIRED)
# Find gslib
find_package(gslib REQUIRED)

# Sources
file(GLOB LIB_SOURCES  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
file(GLOB EXAMPLE_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/example/*.c")
file(GLOB TEST_SRCS    "${CMAKE_CURRENT_SOURCE_DIR}/test/*.c")

# Library target
if(PARRSB_SHARED)
  add_library(parRSB SHARED ${LIB_SOURCES})
else()
  add_library(parRSB STATIC ${LIB_SOURCES})
endif()

set_target_properties(parRSB PROPERTIES OUTPUT_NAME parRSB)

# Include directories
target_include_directories(parRSB PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Compile flags
target_compile_options(parRSB PRIVATE
  -std=c99
  -Wall -Wextra -Wpedantic
  -Wno-unused-function
  -Wno-c23-extensions
)

# Preprocessor definitions
if(PARRSB_SYNC STREQUAL "reduction")
  target_compile_definitions(parRSB PRIVATE PARRSB_SYNC_REDUCTION)
elseif(PARRSB_SYNC STREQUAL "barrier")
  target_compile_definitions(parRSB PRIVATE PARRSB_SYNC_BARRIER)
endif()

if(NOT PARRSB_VISIBILITY STREQUAL "")
  target_compile_definitions(parRSB PRIVATE PARRSB_VISIBILITY=${PARRSB_VISIBILITY})
endif()

# Link libraries
target_link_libraries(parRSB PRIVATE gs::gs MPI::MPI_C m)

install(TARGETS parRSB
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.h"
)

# Shared helper for examples & tests
# Collects common compile/link settings so they don't have to be repeated.
function(add_mpi_executable target_name source_file extra_include)
  add_executable(${target_name} "${source_file}")

  target_include_directories(${target_name} PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${GSLIB_DIR}/include"
    ${extra_include}
  )

  target_compile_options(${target_name} PRIVATE
    -std=c99
    -Wall -Wextra -Wpedantic
    -Wno-unused-function
    -Wno-c23-extensions
    $<IF:$<BOOL:${DEBUG}>,-g,-O2>
  )

  target_link_libraries(${target_name} PRIVATE
    parRSB
    MPI::MPI_C
    "${GSLIB_DIR}/lib/libgs.a"
    m
  )

  add_dependencies(${target_name} parRSB)
endfunction()

# Example binaries
foreach(src ${EXAMPLE_SRCS})
  get_filename_component(name "${src}" NAME_WE)
  add_mpi_executable("example_${name}" "${src}" "")
endforeach()

# Test binaries
enable_testing()

foreach(src ${TEST_SRCS})
  get_filename_component(name "${src}" NAME_WE)
  add_mpi_executable("test_${name}" "${src}"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/test"
  )
  add_test(NAME "${name}" COMMAND "test_${name}")
endforeach()

# Format target (optional, requires clang-format)
find_program(CLANG_FORMAT clang-format)
if(CLANG_FORMAT)
  file(GLOB_RECURSE ALL_C_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/example/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/test/*.c"
  )
  add_custom_target(format
    COMMAND ${CLANG_FORMAT} -i ${ALL_C_FILES}
    COMMENT "Running clang-format"
  )
endif()

# Configure library export options
set(EXPORT_TARGETS_FILE "parRSBTargets.cmake")
# Make it exportable by creating an export name
install(TARGETS parRSB
  EXPORT parRSBTargets
  ARCHIVE DESTINATION lib
)
install(EXPORT parRSBTargets
  FILE parRSBTargets.cmake
  NAMESPACE parRSB::
  DESTINATION lib/cmake/parRSB
)
# Install headers and CMake configuration
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/parRSBConfigVersion.cmake"
  VERSION 1.0.9
  COMPATIBILITY AnyNewerVersion
)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/parRSBConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/parRSBConfig.cmake"
  @ONLY
)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/parRSBConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/parRSBConfigVersion.cmake"
  DESTINATION lib/cmake/parRSB
)
install(EXPORT parRSBTargets
  FILE parRSBTargets.cmake
  NAMESPACE parRSB::
  DESTINATION lib/cmake/parRSB
)

file(GLOB HEADER_FILES "${SRCDIR}/*.h")
install(FILES ${HEADER_FILES} DESTINATION include)
