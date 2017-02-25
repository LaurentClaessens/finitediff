####### Compiler, tools and options

CXX           = LC_ALL=C g++ -std=c++17
DEL_FILE      = rm -f

BUILD_DIR = build/
SRC_DIR = src/

####### Implicit rules

clean: 
	-$(DEL_FILE) *.o
	-$(DEL_FILE) *~ core *.core


####### Compile
all:  finitediff
RepeatFunction : $(SRC_DIR)RepeatFunction.cpp $(SRC_DIR)RepeatFunction.h
	$(CXX) -c  -o $(BUILD_DIR)RepeatFunction.o $(SRC_DIR)RepeatFunction.cpp
finitediff: RepeatFunction
	$(CXX) -g  -o finitediff $(BUILD_DIR)RepeatFunction.o $(SRC_DIR)finitediff.cpp
