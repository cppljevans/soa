#ifndef SSE_HPP_INCLUDED
#define SSE_HPP_INCLUDED
#if defined(__GNUC__)
  #include <emmintrin.h>
#else
  #include <mmintrin.h>
#endif//for __m128 type and _mm_* functions.
#include <cstddef>//for std::size_t
using sse_float=float;
unsigned constexpr sse_width=4;//"width" of vector operands in single sse instruction.
std::size_t constexpr sse_align=sse_width*sizeof(sse_float);//alignment required for sse operand instructions.
#endif// SSE_HPP_INCLUDED
