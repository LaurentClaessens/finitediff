#! /bin/bash


# One script for
#
# - clean 
# - recall your git status
# - compile all
# - tests

mkdir build
git status
make clean

echo
echo "TESTS --------------------"
echo

function launch_test 
{
    echo "+++ Compilation : " $1
    make $1
    ./build/$1
    echo "--- Ended " $1
}

launch_test "sn_permutation_unit_tests"
launch_test "gauss_unit_tests"
launch_test "exceptions_unit_tests"
launch_test "repeat_function_unit_tests"
launch_test "sn_multiplication_unit_tests"
launch_test "multiplication_unit_tests"
launch_test "sn_matrix_unit_tests"
launch_test "sn_line_unit_tests"
launch_test "sn_element_unit_tests"
launch_test "plu_unit_tests"

echo
echo "COMPILATION --------------------"
echo

make all

