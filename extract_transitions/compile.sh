#!/bin/bash

##gcc -Wall -g -o extract_transitions extract_transitions.c knhx.c kstring.c
gcc -Wall -O2 -o extract_transitions extract_transitions.c knhx.c kstring.c
