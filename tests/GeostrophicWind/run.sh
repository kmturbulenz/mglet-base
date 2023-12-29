#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    # running the executable
    mpirun -n 1 $MGLET_BIN 2>&1 | tee mglet.OUT
    # conducting a check of the results (uses "LOGS/uvwbulk.log")
    python check.py
elif [[ "$ACTION" == "clean" ]]; then
    # removing generated files
    rm -rf LOGS fields.h5 mglet-perf-report.txt *.OUT
else
    # response to undefined action
    echo "Invalid action: $ACTION"
    exit 1
fi
