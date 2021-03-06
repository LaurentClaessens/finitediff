#! /bin/bash
# -*- coding: utf8 -*-


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

CURRENT_DIR=`pwd`
DEPLOY_LOG_FILE=$CURRENT_DIR/.deploy.log
TEST_LOG_FILE=$CURRENT_DIR/.testing.log

rm $DEPLOY_LOG_FILE
touch $DEPLOY_LOG_FILE
rm $TEST_LOG_FILE
touch $TEST_LOG_FILE


mkdir build > /dev/null
git status  > $DEPLOY_LOG_FILE
make clean

./testing.sh $TEST_LOG_FILE

echo "Deploy results -------------"
cat $DEPLOY_LOG_FILE
