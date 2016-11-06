#ifndef BIT_VECTOR_HPP_INCLUDED
#define BIT_VECTOR_HPP_INCLUDED
#include <cstdint>//uint64_t
#include <cstddef>//std::size_t
#include "bitsof.hpp"
//#define HAVE_GOON_BIT_VECTOR
#ifdef HAVE_GOON_BIT_VECTOR
  std::size_t constexpr bits_per_block
    =bitsof<uint64_t>()//must be multiple of sse_width.
    ;  
  #include <goon/bit_vector.hpp>
  using goon::bit_vector;
#else
  std::size_t constexpr bits_per_block
    =bitsof<uint64_t>()//must be multiple of sse_width.
    ;  
  #include <vector>
    //The following bit_vector implementation
    //was copy&pasted from:
    /*
  From: Michael Marcin <mike.marcin@gmail.com>
  Newsgroups: gmane.comp.lib.boost.devel
  Subject: Re: interest in structure of arrays container?
  Date: Sun, 30 Oct 2016 18:45:41 -0500
     */
    //and then minor modifications made.
    //Then converted to use bitset instead of uint64_t
    //as value_type for convenience.
  #include <bitset>
  struct bit_vector 
  : public std::vector
    < std::bitset<bits_per_block>
    >
  {
      using 
    super_t=
      std::vector
      < std::bitset<bits_per_block>
      >;
    using super_t::data;
    
    bit_vector():num_bits_m(0){}  
    void resize( std::size_t num_bits, bool value = false )
    {
        std::size_t num_blocks = 
          ( num_bits + (bits_per_block-1) ) 
          / bits_per_block
          ;//assert(num_blocks*bits_per_block >= num_bits);
        super_t::value_type block;//all bits false.
        if(value) block.set();//set all bits true.
        super_t::resize( num_blocks, block );
        num_bits_m = num_bits;
    }
    std::size_t size() const { return num_bits_m; }
    std::size_t vec_size() const { return super_t::size(); }
    auto* end(){ return data()+vec_size();}
  private:  
    std::size_t num_bits_m;
  };   
#endif//HAVE_GOON_BIT_VECTOR  
#endif//BIT_VECTOR_HPP_INCLUDED
