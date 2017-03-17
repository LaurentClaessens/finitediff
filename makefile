####### Compiler, tools and options

CXX           = LC_ALL=C g++ -std=c++17
CXXFLAGS      = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)


DEL_FILE      = rm -f

BUILD_DIR = build/
SRC_DIR = src/
TESTS_DIR = tests/




# NOTE :
# 
# SNmatrix being a template class, it does not need to be explicitly compiled here.
#     More : trying to link it in 'unit_tests' causes the error
#     SNmatrix.o: file format not recognized; treating as linker script 

####### Compile
all:  finitediff unit_tests RepeatFunction 
clean:
	$(DEL_FILE) build/*
RepeatFunction : $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(CXX) $(CXXFLAGS)  -c   $(SRC_DIR)RepeatFunction.cpp   -o $(BUILD_DIR)RepeatFunction.o
finitediff: RepeatFunction $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(CXX) $(CXXFLAGS)   $(SRC_DIR)finitediff.cpp $(BUILD_DIR)RepeatFunction.o  -o $(BUILD_DIR)finitediff

exceptions_unit_tests: $(TESTS_DIR)exceptions_unit_tests.cpp
	$(CXX) $(CXXFLAGS)   $(TESTS_DIR)exceptions_unit_tests.cpp    -lcppunit -o $(BUILD_DIR)exceptions_unit_tests

multiplication_unit_tests: $(TESTS_DIR)multiplication_unit_tests.cpp
	$(CXX) $(CXXFLAGS)   $(TESTS_DIR)multiplication_unit_tests.cpp    -lcppunit -o $(BUILD_DIR)multiplication_unit_tests
	
sn_matrix_tests: $(TESTS_DIR)sn_matrix_unit_tests.cpp
	$(CXX) $(CXXFLAGS)   $(TESTS_DIR)sn_matrix_unit_tests.cpp    -lcppunit -o $(BUILD_DIR)sn_matrix_unit_tests

repeat_function_unit_tests: RepeatFunction  $(TESTS_DIR)repeat_function_unit_tests.cpp
	$(CXX) $(CXXFLAGS)   $(TESTS_DIR)repeat_function_unit_tests.cpp $(BUILD_DIR)RepeatFunction.o  -lcppunit -o $(BUILD_DIR)repeat_function_unit_tests
	
unit_tests: repeat_function_unit_tests exceptions_unit_tests multiplication_unit_tests  sn_matrix_tests
