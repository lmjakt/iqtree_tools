#!/bin/bash

if [[ $1 == "debug" ]]
then
    ## for testing:
    echo "gcc -Wall -g -o extract_transitions extract_transitions.c knhx.c kstring.c"
    gcc -Wall -g -o extract_transitions extract_transitions.c knhx.c kstring.c
else
    echo "gcc -Wall -O2 -o extract_transitions extract_transitions.c knhx.c kstring.c"
    gcc -Wall -O2 -o extract_transitions extract_transitions.c knhx.c kstring.c
fi
