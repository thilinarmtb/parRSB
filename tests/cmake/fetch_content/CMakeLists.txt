cmake_minimum_required(VERSION 3.24)

project(parRSB_test LANGUAGES C)

include(FetchContent)
FetchContent_Declare(
  parRSB
  SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../"
  FIND_PACKAGE_ARGS
)
FetchContent_MakeAvailable(parRSB)

add_executable(exe "${CMAKE_CURRENT_SOURCE_DIR}/src.c")
target_link_libraries(exe PRIVATE parRSB)
install(TARGETS exe)
