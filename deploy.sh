#! /bin/bash


# One script for
#
# - clean 
# - recall your git status
# - compile all
# - tests

function launch_test 
{
    echo "** Launching" $1
    ./build/$1
}

mkdir build
git status
make clean

echo
echo "COMPILATION --------------------"
echo

make all

echo
echo "TESTS --------------------"
echo

launch_test "exceptions_unit_tests"
launch_test "repeat_function_unit_tests"
launch_test "multiplication_unit_tests"
