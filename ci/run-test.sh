#!/bin/bash

INPUTDIR=$1
MGLET_BIN=$2

# This will print a backtrace if Ctrl-C (SIGINT) is received
# Might not work as expected if several processes with the same binary are
# running at the same computer at the same time (will give the backtrace
# of all processes with the same name). Since we cannot attach to the PID
# of the mpirun command (which is different from the PID of each
# individual MGLET process) this is a good shot.
print_backtrace()
{
    for PID in $(pidof $MGLET_BIN); do
        echo "Backtrace of $MGLET_BIN process: $PID:"
        gdb -batch -ex "attach $PID" -ex bt
        echo ""
    done
    exit 1
}
trap print_backtrace SIGINT

# Copy input files to present directory
cp $INPUTDIR/* .

# Run testcase
TIC=`date +%s`
./run.sh test $MGLET_BIN &
wait $!
RETCODE=$?
TOC=`date +%s`
RUNTIME=$((TOC-TIC))

if [[ -v GITHUB_STEP_SUMMARY ]]; then
    if [[ ! -s $GITHUB_STEP_SUMMARY ]]; then
        echo "# Test summary" >> $GITHUB_STEP_SUMMARY
        echo "" >> $GITHUB_STEP_SUMMARY
        echo "| Case | Status | Duration (s) |" >> $GITHUB_STEP_SUMMARY
        echo "| ---- | ------ | -----------  |" >> $GITHUB_STEP_SUMMARY
    fi

    CASENAME=${PWD##*/}

    if [[ $RETCODE -eq "0" ]]; then
        RESULT=":green_circle: Success"
    else
        RESULT=":red_circle: Failure"
    fi

    echo "| $CASENAME | $RESULT | $RUNTIME |" >> $GITHUB_STEP_SUMMARY
fi

exit $RETCODE
