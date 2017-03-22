####### Compiler, tools and options

CXX           = LC_ALL=C g++ -std=c++17
CLANG         = LC_ALL=C clang++ -std=c++14

COMPILATOR = $(CLANG)

CXXFLAGS      = -pipe -O2 -Wall -W -D_REENTRANT $(DEFINES)


DEL_FILE      = rm -f

BUILD_DIR = build/
SRC_DIR = src/
SNMATRICES_DIR = $(SRC_DIR)/SNmatrices/
TESTS_DIR = tests/


# NOTE for myself (because I'm a noob) :
# 
# SNmatrix being a template class, it does not need to be explicitly compiled here.
#     More : trying to link it in 'unit_tests' causes the error
#     SNmatrix.o: file format not recognized; treating as linker script 

####### Compile
all:  finitediff unit_tests RepeatFunction 
clean:
	$(DEL_FILE) build/*
RepeatFunction : $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(COMPILATOR) $(CXXFLAGS)  -c   $(SRC_DIR)RepeatFunction.cpp   -o $(BUILD_DIR)RepeatFunction.o
finitediff: RepeatFunction $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(COMPILATOR) $(CXXFLAGS)   $(SRC_DIR)finitediff.cpp $(BUILD_DIR)RepeatFunction.o  -o $(BUILD_DIR)finitediff
repeat_function_unit_tests: RepeatFunction  $(TESTS_DIR)repeat_function_unit_tests.cpp
	$(COMPILATOR) $(CXXFLAGS)   $(TESTS_DIR)$@.cpp $(BUILD_DIR)RepeatFunction.o  -lcppunit -o $(BUILD_DIR)$@

m_num: $(SNMATRICES_DIR)m_num.cpp  $(SNMATRICES_DIR)m_num.h
	$(COMPILATOR) $(CXXFLAGS)  -c   $(SNMATRICES_DIR)$@.cpp   -o $(BUILD_DIR)$@.o

exceptions_unit_tests: $(TESTS_DIR)exceptions_unit_tests.cpp m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp   $(BUILD_DIR)m_num.o     -lcppunit -o $(BUILD_DIR)$@

multiplication_unit_tests: $(TESTS_DIR)multiplication_unit_tests.cpp m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp   $(BUILD_DIR)m_num.o     -lcppunit -o $(BUILD_DIR)$@
	
sn_matrix_unit_tests: $(TESTS_DIR)sn_matrix_unit_tests.cpp m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp   $(BUILD_DIR)m_num.o     -lcppunit -o $(BUILD_DIR)$@

sn_line_unit_tests: $(TESTS_DIR)sn_line_unit_tests.cpp m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp   $(BUILD_DIR)m_num.o     -lcppunit -o $(BUILD_DIR)$@

sn_element_unit_tests: $(TESTS_DIR)sn_element_unit_tests.cpp m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp   $(BUILD_DIR)m_num.o     -lcppunit -o $(BUILD_DIR)$@

gauss_unit_tests: $(TESTS_DIR)gauss_unit_tests.cpp  m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp   $(BUILD_DIR)m_num.o     -lcppunit -o $(BUILD_DIR)$@

plu_unit_tests: $(TESTS_DIR)plu_unit_tests.cpp  m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp   $(BUILD_DIR)m_num.o     -lcppunit -o $(BUILD_DIR)$@

sn_permutation_unit_tests: $(TESTS_DIR)sn_permutation_unit_tests.cpp  m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp  $(BUILD_DIR)m_num.o  -lcppunit -o $(BUILD_DIR)$@

sn_multiplication_unit_tests: $(TESTS_DIR)sn_permutation_unit_tests.cpp  m_num
	$(COMPILATOR) $(CXXFLAGS) -g  $(TESTS_DIR)$@.cpp  $(BUILD_DIR)m_num.o  -lcppunit -o $(BUILD_DIR)$@
	
unit_tests: m_num repeat_function_unit_tests exceptions_unit_tests multiplication_unit_tests sn_matrix_unit_tests\
	sn_line_unit_tests sn_element_unit_tests gauss_unit_tests plu_unit_tests sn_multiplication_unit_tests sn_permutation_unit_tests
