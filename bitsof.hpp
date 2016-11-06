#ifndef BITSOF_HPP_INCLUDED
#define BITSOF_HPP_INCLUDED
//Purpose:
//  return bits in a type.
//Reference:
/*
http://www.cplusplus.com/reference/climits/
 */
//================ 
#include <climits>//CHAR_BIT
template<typename T>
constexpr auto bitsof(){ return CHAR_BIT*sizeof(T);}
#endif//BITSOF_HPP_INCLUDED
