#! /bin/bash
# -*- coding: utf8 -*-

TEST_LOG_FILE=$1

echo
echo "TESTS --------------------"
echo

echo $TEST_LOG_FILE
exit

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
            echo -e "---- ${RED} The test ${CYAN}" $1 "${RED} got a problem ${NC}" >> $TEST_LOG_FILE
            echo "----"
    fi
    echo "--- Ended " $1
}

launch_test "include_plu_tests"
launch_test "plu_unit_tests"
launch_test "multigauss_unit_tests"
launch_test "m_num_unit_tests"
launch_test "sn_permutation_unit_tests"
launch_test "exceptions_unit_tests"
launch_test "repeat_function_unit_tests"
launch_test "sn_multiplication_unit_tests"
launch_test "multiplication_unit_tests"
launch_test "sn_matrix_unit_tests"
launch_test "sn_line_unit_tests"
launch_test "sn_element_unit_tests"
launch_test "gauss_unit_tests"
launch_test "utilities_tests"
launch_test "sn_gaussian_unit_tests"

echo "Tests results ---------------"
echo "In $TEST_LOG_FILE"
cat $TEST_LOG_FILE

