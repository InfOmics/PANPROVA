#!/bin/bash

# example use:
# bash compile.sh <optional compiler flags>
# bash compile.sh -DVERBOSE : enable verbose logging
# bash compile.sh -DDEBUG
# bash compile.sh -DVERBOSE -DDEBUG : enable both verbose logging and debug mode

g++ create_hgt_pool.cpp -o create_hgt_pool "$@"
g++ evolve.cpp -o evolve "$@"
