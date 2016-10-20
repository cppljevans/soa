#ifndef VEC_OFFSETS_HPP_INCLUDED
#define VEC_OFFSETS_HPP_INCLUDED
#include <cassert>
#include <array>
#include <cstddef> //std::size_t
//#define TRACE_NXT_OFFSET
#ifdef TRACE_NXT_OFFSET
  #include <iostream>
#endif
  std::size_t
nxt_offset
  ( std::size_t now_offset
  , std::size_t now_size
  , std::size_t nxt_align
  )
  /**@brief
   *  Find the next offset which provides
   *  enough space(now_size) from now_offset
   *  and has alignment, nxt_align.
   */
  {
    std::size_t result=now_offset+now_size;
  #ifdef TRACE_NXT_OFFSET
    std::cout<<__func__
      <<":now_offset="<<now_offset
      <<":now_size="<<now_size
      <<":nxt_align="<<nxt_align
      <<":result1="<<result
      <<"\n";
  #endif
    std::size_t remainder=result%nxt_align;
    if(remainder)
      result+=nxt_align-remainder; //if result not aligned, make it aligned.
  #ifdef TRACE_NXT_OFFSET
    std::cout<<__func__<<":remainder1="<<remainder<<":result2="<<result<<"\n";
  #endif
    remainder=result%nxt_align;
  #ifdef TRACE_NXT_OFFSET
    std::cout<<__func__<<":remainder2="<<remainder<<"\n";
  #endif
    assert(remainder == 0);
    return result;
  }
  template<typename... Ts>
  auto
vec_offsets
  ( std::size_t vec_size
  )
  /**@brief
   *  Calculate offsets in a char buffer for storing 
   *  suitably aligned Ts[vec_size]...
   *  Last offset(at sizeTs) is for the char buffer size.
   *  ***CAUTION*** the last offset should not be used
   *  if 2 or more of such character buffers are contiguous.
   *  That's because the 2nd such buffer may not have the
   *  correct alignment for one or more of the Ts[vec_size]...
   *  The test driver, vec_offsets.test.cpp, shows this.
   */
  {
    std::size_t const sizeTs=sizeof...(Ts);
    std::array<std::size_t,sizeTs+1> result{0};
    std::size_t sizes[sizeTs]={sizeof(Ts)...};
    std::size_t aligns[sizeTs]={alignof(Ts)...};
    for(std::size_t i_T=1; i_T<sizeTs; ++i_T)
    {
      result[i_T]=nxt_offset( result[i_T-1], sizes[i_T-1]*vec_size, aligns[i_T]);
    }
    result[sizeTs]=result[sizeTs-1]+sizes[sizeTs-1]*vec_size;
    //^The size of buffer to hold all the Ts[vec_size]....
    return result;
  }
#endif//VEC_OFFSETS_HPP_INCLUDED
