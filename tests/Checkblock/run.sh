#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    mpirun -n 2 $MGLET_BIN parameters-checkblock.json 2>&1 | tee checkblock.OUT
    mpirun -n 2 $MGLET_BIN parameters-newblock.json 2>&1 | tee newblock.OUT
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS grids_new.h5 ib_stencils.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action: $ACTION"
    exit 1
fi
