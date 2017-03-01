####### Compiler, tools and options

CXX           = LC_ALL=C g++ -std=c++17
DEL_FILE      = rm -f

CPPUNIT_SO = /usr/lib/i386-linux-gnu/libcppunit.so
CPPUNIT_A = /usr/lib/i386-linux-gnu/libcppunit.a

BUILD_DIR = build/
SRC_DIR = src/
TESTS_DIR = tests/


####### Compile
all:  finitediff tests
clean:
	$(DEL_FILE) build/*.o
	$(DEL_FILE) tests_RepeatFunction
	$(DEL_FILE) finitediff
RepeatFunction : $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(CXX) -c  -o $(BUILD_DIR)RepeatFunction.o $(SRC_DIR)RepeatFunction.cpp
finitediff: RepeatFunction $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(CXX) -g  $(SRC_DIR)finitediff.cpp $(BUILD_DIR)RepeatFunction.o  -o finitediff
tests: RepeatFunction  $(TESTS_DIR)test_RepeatFunction.cpp 
	$(CXX) -g $(TESTS_DIR)test_RepeatFunction.cpp  $(BUILD_DIR)RepeatFunction.o -lcppunit   -o tests_RepeatFunction
