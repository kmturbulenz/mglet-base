#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    $MGLET_BIN parameters-ramp.json 2>&1 | tee mglet-ramp.OUT
    $MGLET_BIN parameters-ramp_inf.json 2>&1 | tee mglet-ramp_inf.OUT

    if python3 &> /dev/null
    then
        python3 plot.py 2>&1 | tee plot.OUT
    fi
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields.h5 ib_stencils.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action: $ACTION"
    exit 1
fi
