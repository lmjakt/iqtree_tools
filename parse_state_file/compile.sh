#!/bin/bash

if [[ $1 == "debug" ]]
then
    ## for testing:
    echo "gcc -Wall -g -o state_to_nibble state_to_nibble.c kstring.c -lz"
    gcc -Wall -g -o state_to_nibble state_to_nibble.c kstring.c -lz
else
    ## for release
    echo "gcc -Wall -O2 -o state_to_nibble state_to_nibble.c kstring.c -lz"
    gcc -Wall -O2 -o state_to_nibble state_to_nibble.c kstring.c -lz
fi

