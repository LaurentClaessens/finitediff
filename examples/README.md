## Examples

This directory contains examples of the functionaries developed so far. This is only high level functionalities.

* PLU decomposition.

### Usage

Before to launch the examples, you should try to perform the tests, to be sure that everything if okay. For that :
```
./deploy.sh
```
Now for compiling the examples :
```
make examples
```

Note : you have to create the subdirectory `build` (already there if you launched the tests).

The executable is created in the subdirectory `build.`

### In your own program

#### PLU in you program

The main have to include the following :
```
#include "../src/SNplu.h"
```

For compiling : 
```
clang++ -std=c++14 -pipe -O2 -Wall -W -D_REENTRANT   -g  YOUR_SOURCE_CPP_FILE   build/m_num.o  -o YOUR_TARGET_BUILD_FILE
```

