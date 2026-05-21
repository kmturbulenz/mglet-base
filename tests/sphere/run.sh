#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    mpirun -n 1 $MGLET_BIN parameters-1.json 2>&1 | tee mglet.OUT
    mpirun -n 1 $MGLET_BIN parameters-2.json 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields-checkpoint.h5 fields.h5 ib_stencils.h5 probes.h5 snapshots.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action: $ACTION"
    exit 1
fi
