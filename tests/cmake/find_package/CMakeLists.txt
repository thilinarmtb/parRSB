cmake_minimum_required(VERSION 3.24)

project(parRSB_test LANGUAGES C)

find_package(parRSB REQUIRED)

add_executable(exe "${CMAKE_CURRENT_SOURCE_DIR}/src.c")
target_link_libraries(exe PRIVATE parRSB::parRSB)
install(TARGETS exe)
