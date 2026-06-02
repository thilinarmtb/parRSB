#!/bin/bash
set -e
set -o pipefail

: ${CC:=mpicc}
: ${CMAKE:=cmake}
: ${MPIRUN:=mpirun}
: ${MPIOPTS:=}
: ${NP:="1 2 3 4"}
: ${CLANG_FORMAT:=clang-format}

############################
# Don't touch what follows #
############################
WRK_DIR=$(mktemp -d)
TST_DIR=${WRK_DIR}/tests
BLD_DIR=${WRK_DIR}/build
INS_DIR=${WRK_DIR}/install
NEK_DIR=${WRK_DIR}/nek5k_ci
GS_SRC_DIR=${WRK_DIR}/gslib
GS_BLD_DIR=${GS_SRC_DIR}/build
GS_INS_DIR=${GS_SRC_DIR}/install

SRC_FILES=

###################
# Redirect output #
###################
exec 4>&1
exec  >${WRK_DIR}/out.log

####################
# Helper functions #
####################
function msg() {
  printf "%s" "$1" >&4
}

function msgln() {
  printf "%s\n" "$1" >&4
}

########################################
# Functions to help build gslib/parRSB #
########################################
function build_gslib() {
  msg "Building gslib ... "
  git clone -q https://github.com/thilinarmtb/gslib.git ${GS_SRC_DIR}
  cmake -B ${GS_BLD_DIR} -S ${GS_SRC_DIR} -DCMAKE_INSTALL_PREFIX=${GS_INS_DIR}
  cmake --build ${GS_BLD_DIR} --target install
  msgln "ok"
}

function build_parrsb() {
  msg "Building parRSB ... "
  CC=${CC} ${CMAKE} -B ${BLD_DIR} -S . -DCMAKE_INSTALL_PREFIX=${INS_DIR} \
    -Dgs_DIR=${GS_INS_DIR}/lib/cmake/gs
  ${CMAKE} --build ${BLD_DIR} --target install
  msgln "ok"
}

function setup_format() {
  if ! command -v clang-format >/dev/null 2>&1; then
    msgln "clang-format not found in PATH!"
    exit 1
  fi

  SCRIPT_DIR=$(dirname -- "$(realpath -- "${BASH_SOURCE[0]}")")
  SRC_DIR=${SCRIPT_DIR}/../
  SRC_FILES=$(find ${SRC_DIR} -type f -name "*.[ch]")
}

function format() {
  msg "Running clang-format ... "
  setup_format
  ${CLANG_FORMAT} -i ${SRC_FILES}
  msgln "ok"
}

#############
# Run tests #
#############
function check_format() {
  msg "Running clang-format checks ... "
  setup_format
  ${CLANG_FORMAT} --dry-run -Werror -i ${SRC_FILES}
  msgln "ok"
}

function test_cmake_exported_target() {
  cd ${TST_DIR}/cmake/find_package
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX=./install \
    -DparRSB_DIR=${INS_DIR}/lib/cmake/parRSB
  cmake --build build --target install
}

function test_cmake_inclusion_with_fetch_content() {
  cd ${TST_DIR}/cmake/find_package
  cmake -S . -B build -DCMAKE_INSTALL_PREFIX=./install
  cmake --build build --target install
}

function cmake_tests() {
  msg "Running cmake tests ... "
  mkdir -p ${TST_DIR}
  cp -r tests/* ${TST_DIR}
  test_cmake_exported_target
  test_cmake_inclusion_with_fetch_content
  msgln "ok"
}

function unit_tests() {
  msgln "Running unit tests ..."
  for test_bin in `ls ${INS_DIR}/bin/[0-9][0-9][0-9]_*`; do
    for np in ${NP}; do
      msg "  NP=${np}, $(basename ${test_bin}) "
      ${MPIRUN} ${MPIOPTS} -np ${np} ${test_bin}
      msgln "... ok"
    done
  done
}

function nek5k_tests() {
  meshes=(box_2x2x2 pyramid tgv e3q solid ethier vortex expansion)

  git clone -q https://github.com/thilinarmtb/parRSB-github-ci.git ${NEK_DIR}

  export PARRSB_RSB_ALGO=0
  export PARRSB_VERBOSE_LEVEL=2

  msgln "Running gencon tests ..."
  for mesh in "${meshes[@]}"; do
    cd ${NEK_DIR}/${mesh}
    tol=(`cat test.txt | grep tol`); tol=${tol[2]}

    for np in ${NP}; do
      msg "  NP=${np}, ${mesh} "
      ${MPIRUN} ${MPIOPTS} -np ${np} ${INS_DIR}/bin/gencon --mesh ${mesh} \
            --tol=${tol} --dump=0 --test=1
      msgln "... ok"
    done
  done

  err_cnt=0
  msgln "Running genmap tests ..."
  # there are some allowed failures
  set +e
  for mesh in "${meshes[@]}"; do
    cd ${NEK_DIR}/${mesh}
    tol=(`cat test.txt | grep tol`); tol=${tol[2]}

    for np in ${NP}; do
      msg "  NP=${np}, ${mesh} "
      ${MPIRUN} ${MPIOPTS} -np ${np} ${INS_DIR}/bin/genmap --mesh ${mesh} \
        --tol=${tol} --dump=0 --test=1
      if [[ $? -ne 0 ]]; then
        err_cnt=$(( err_cnt + 1 ))
        msgln "... not ok"
      else
        msgln "... ok"
      fi
    done
  done
  set -e

  if [[ $err_cnt -gt 5 ]]; then
    exit 1
  fi
}

if [[ "$1" == "-fmt" ]]; then
  format
  exit 0
elif [[ "$1" == "-check-fmt" ]]; then
  check_format
  exit 0
fi
msgln "Running tests in ${WRK_DIR} ..."
msgln "  CC    : ${CC}"
msgln "  cmake : `which ${CMAKE}`"
msgln "  mpirun: `which ${MPIRUN}`, opts: \"${MPIOPTS}\", np: ${NP}"
build_gslib
build_parrsb
cmake_tests
unit_tests
nek5k_tests
msgln "Tests passed."
