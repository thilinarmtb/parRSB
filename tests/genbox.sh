#!/bin/bash

NEK5K=${1:-"${HOME}/Nek5000"}
MESH=${2:-"./linear_periodic_mesh.box"}

GENBOX=${NEK5K}/bin/genbox
if [[ ! -f $GENBOX ]]; then
    echo "Error: genbox not found in \"$GENBOX\""
    exit 1
fi

abs_path=$(realpath $MESH)
$GENBOX << EOF
${abs_path}
EOF
