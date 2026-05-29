#!/bin/bash
set -e
set -o pipefail


: ${CC:=mpicc}
: ${CMAKE:=cmake}
: ${MPIRUN:=mpirun}
: ${MPIOPTS:=}
: ${NP:=4}

############################
# Don't touch what follows #
############################
WRK_DIR=$(mktemp -d)
TST_DIR=${WRK_DIR}/tests
BLD_DIR=${WRK_DIR}/build
INS_DIR=${WRK_DIR}/install
GS_SRC_DIR=${WRK_DIR}/gslib
GS_BLD_DIR=${GS_SRC_DIR}/build
GS_INS_DIR=${GS_SRC_DIR}/install

########################################
# Functions to help build gslib/parRSB #
########################################
function build_gslib() {
  git clone https://github.com/thilinarmtb/gslib.git ${GS_SRC_DIR}
  cmake -B ${GS_BLD_DIR} -S ${GS_SRC_DIR} -DCMAKE_INSTALL_PREFIX=${GS_INS_DIR}
  cmake --build ${GS_BLD_DIR} --target install
}

function build_parrsb() {
  CC=${CC} ${CMAKE} -B ${BLD_DIR} -S . -DCMAKE_INSTALL_PREFIX=${INS_DIR} \
    -Dgs_DIR=${GS_INS_DIR}/lib/cmake/gs
  ${CMAKE} --build ${BLD_DIR} --target install
}

#############
# Run tests #
#############
function check_formatting() {
  echo "Running format checks ..."
  ${CMAKE} --build ${BLD_DIR} --target format-check
}

function test_cmake_exported_target() {
  cd ${TST_DIR}/cmake/find_package
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX=${INS_DIR} \
    -DparRSB_DIR=${INS_DIR}/lib/cmake/parRSB
  cmake --build build --target install
  cd -
}

function test_cmake_inclusion_with_fetch_content() {
  cd ${TST_DIR}/cmake/find_package
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX=${INS_DIR}
  cmake --build build --target install
  cd -
}

function run_cmake_tests() {
  echo "Running cmake tests ..."
  mkdir -p ${TST_DIR}
  cp -r tests/* ${TST_DIR}
  test_cmake_exported_target
  test_cmake_inclusion_with_fetch_content
}

function run_unit_tests() {
  echo "Running unit tests ..."
  for test_bin in `ls ${INS_DIR}/bin/[0-9][0-9][0-9]_*`; do
    echo -n "Running test: ${test_bin}"
    ${MPIRUN} ${MPIOPTS} -np ${NP} -- ${test_bin}
    echo " ... ok"
  done
}

GREEN='\033[0;32m'
RESET='\033[0m'

echo -e "${GREEN}"
echo "Running tests in ${WRK_DIR} ..."
echo "  CC    : ${CC}"
echo "  cmake : `which ${CMAKE}`"
echo "  mpirun: `which ${MPIRUN}`, opts: \"${MPIOPTS}\", np: ${NP}"
echo -e "${RESET}"
build_gslib
build_parrsb
echo -e "${GREEN}"
check_formatting
run_cmake_tests
run_unit_tests
echo -e "Tests passed."
echo -e "${RESET}"
