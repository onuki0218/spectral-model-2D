#!/bin/bash

rm -f sub.sh.*
rm -f sub_restart.sh.*
rm -f debug.sh.*
rm -f debug_restart.sh.*
if [ $# -ne 0 ]; then
    rm -r -f ../data/data$1
fi
