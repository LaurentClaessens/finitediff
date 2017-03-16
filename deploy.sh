#! /bin/bash


# One script for
#
# - clean 
# - compile all
# - test

mkdir build
git status
make clean
make all
./build/exceptions_unit_tests
./build/repeat_function_unit_tests
./build/multiplication_unit_tests
