#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    if python3 &> /dev/null
    then
        python3 run-all.py $MGLET_BIN 2>&1 | tee run-all.OUT
    fi
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields.h5 ib_stencils.h5 mglet-perf-report.txt *.OUT probes-*.h5
else
    echo "Invalid action: $ACTION"
    exit 1
fi
