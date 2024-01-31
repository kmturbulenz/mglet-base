#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    mpirun -n 4 $MGLET_BIN parameters-blocking.json 2>&1 | tee mglet-blocking.OUT
    mpirun -n 4 $MGLET_BIN parameters-run.json 2>&1 | tee mglet-run.OUT
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields.h5 ib_stencils.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action: $ACTION"
    exit 1
fi
