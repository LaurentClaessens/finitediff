
#ifndef __DEBUGPRINT_H_160044_
#define __DEBUGPRINT_H_160044_

#include <iostream>

/**
 So I can put `debug_print<<blahblah` everywhere in my code (for debugging purpose) 
 and found them back with a simple `ack`.
 */ 

std::ostream& debug_print(std::cout);

#endif
