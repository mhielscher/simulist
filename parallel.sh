#!/bin/bash

jobs="+1"
if [[ "$1" == "-j" ]]; then
    shift
    jobs="$1"
    shift
fi
nice -n 19 parallel --progress -j $jobs python simulist.py $1 ::: `seq 1 $2` | sort | uniq -c
