cmake_minimum_required(VERSION 3.21...4.0.0)

project(parRSB_test LANGUAGES C)

include(FetchContent)
FetchContent_Declare(
  parrsb_content
  SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../"
)
FetchContent_MakeAvailable(parrsb_content)

add_executable(exe "${CMAKE_CURRENT_SOURCE_DIR}/src.c")
target_link_libraries(exe PRIVATE parRSB)
install(TARGETS exe)
