#!/bin/bash

MESH=${1:-"./linear_periodic_mesh.box"}
GENBOX=${2:-"$HOME/Nek5000/bin/genbox"}

if [[ ! -f $GENBOX ]]; then
    echo "Error: genbox not found in \"$GENBOX\""
    exit 1
fi

abs_path=$(realpath $MESH)
$GENBOX << EOF
${abs_path}
EOF
