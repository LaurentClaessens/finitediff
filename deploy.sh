#! /bin/bash


# One script for
#
# - clean 
# - compile all
# - test

git status
make clean
make all
./build/unit_tests
