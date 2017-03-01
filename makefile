####### Compiler, tools and options

CXX           = LC_ALL=C g++ -std=c++17
DEL_FILE      = rm -f

BUILD_DIR = build/
SRC_DIR = src/
TESTS_DIR = tests/

CXXFLAGS      = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)


####### Compile
all:  finitediff tests
clean:
	$(DEL_FILE) build/*.o
	$(DEL_FILE) tests_RepeatFunction
	$(DEL_FILE) finitediff
RepeatFunction.o : $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(CXX) -c  -o $(BUILD_DIR)RepeatFunction.o $(SRC_DIR)RepeatFunction.cpp
SNmatrix.o: $(SRC_DIR)SNmatrix.h
	$(CXX) -c  -o $(BUILD_DIR)SNmatrix.o $(SRC_DIR)SNmatrix.h
finitediff: RepeatFunction.o $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(CXX) -g  $(SRC_DIR)finitediff.cpp $(BUILD_DIR)RepeatFunction.o  -o finitediff
tests: test_RepeatFunction test_SNmatrix
test_RepeatFunction :RepeatFunction.o $(TESTS_DIR)test_RepeatFunction.cpp 
	$(CXX) -g $(TESTS_DIR)test_RepeatFunction.cpp  $(BUILD_DIR)RepeatFunction.o -lcppunit   -o tests_RepeatFunction
test_SNmatrix: SNmatrix.o $(TESTS_DIR)test_SNmatrix.cpp  
	$(CXX) -g $(TESTS_DIR)test_SNmatrix.cpp  -lcppunit   -o tests_SNmatrix
