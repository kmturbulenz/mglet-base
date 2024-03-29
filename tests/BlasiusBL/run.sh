#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2

    python3 add-expression.py
    mpirun -n 2 $MGLET_BIN 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action: $ACTION"
    exit 1
fi
