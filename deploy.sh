#! /bin/bash


# One script for
#
# - clean 
# - recall your git status
# - compile all
# - tests

# from http://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux
RED='\033[0;31m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

rm deploy.log
touch deploy.log

mkdir build > /dev/null
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
    if [ $? -eq 0 ]; then
            echo "OK"
    else
            echo "----" 
            echo -e "---- ${RED} The test ${CYAN}" $1 "${RED} got a problem ${NC}"
            echo -e "---- ${RED} The test ${CYAN}" $1 "${RED} got a problem ${NC}" >> deploy.log
            echo "----"
    fi
    echo "--- Ended " $1
}

launch_test "multigauss_unit_tests"
launch_test "plu_unit_tests"
launch_test "sn_gaussian_unit_tests"
launch_test "sn_permutation_unit_tests"
launch_test "exceptions_unit_tests"
launch_test "repeat_function_unit_tests"
launch_test "sn_multiplication_unit_tests"
launch_test "multiplication_unit_tests"
launch_test "sn_matrix_unit_tests"
launch_test "sn_line_unit_tests"
launch_test "sn_element_unit_tests"
launch_test "gauss_unit_tests"

cat deploy.log
