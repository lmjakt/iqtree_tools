#!/bin/bash

## for testing:
##gcc -Wall -g -o state_to_nibble state_to_nibble.c kstring.c -lz

## for release
gcc -Wall -O2 -o state_to_nibble state_to_nibble.c kstring.c -lz
