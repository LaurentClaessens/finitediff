
#ifndef __DEBUGPRINT_H_160044_
#define __DEBUGPRINT_H_160044_

#include <iostream>

/**
 So I can put `debug_print<<blahblah` everywhere in my code (for debugging purpose) 
 and found them back with a simple `ack`.
 */ 

std::ostream& debug_print(std::cout);

template <class M>
void debug_matrix_print(const M& mtr,const std::string& name)
{
    std::cout<<name<<std::endl;
    std::cout<<mtr<<std::endl;
}

#endif
