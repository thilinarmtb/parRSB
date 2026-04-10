#!/bin/bash
set -e
set -o pipefail

: ${CC:=mpicc}
: ${CMAKE:=cmake}

############################
# Don't touch what follows #
############################
WORKDIR=$(mktemp -d)
INSTDIR=${WORKDIR}/parRSB/install

#############################################
# Functions to help with parRSB tests       #
#############################################
function build_parrsb() {
  cmake -S . -B ${WORKDIR}/parRSB/build -DCMAKE_INSTALL_PREFIX=${INSTDIR}
}

function test_cmake_find_package() {
  CC=${CC} ${CMAKE} -S ./test/cmake/find_package -B ${WORKDIR}/find_package/build \
    -DCMAKE_INSTALL_PREFIX=${WORKDIR}/find_package/install \
    -DparRSB_DIR=${INSTDIR}/lib/cmake/parRSB
  cmake --build build --target install
}

function test_cmake_fetch_content() {
  ${CC} ${CMAKE} -S ./test/cmake/fetch_content -B ${WORKDIR}/fetch_content/build \
    -DCMAKE_INSTALL_PREFIX=${WORKDIR}/fetch_content/install
  cmake --build build --target install
}

function test_cmake() {
  test_cmake_find_package
  test_cmake_fetch_content
}

GREEN='\033[0;32m'
RESET='\033[0m'

echo -e "${GREEN}"
echo "Running tests in ${WORKDIR} ..."
echo "  CC    : ${CC}"
echo "  cmake : `which ${CMAKE}`"
echo -e "${RESET}"
build_parrsb
test_cmake
echo -e "${GREEN}Tests passed.${RESET}"
