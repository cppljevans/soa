#ifndef BIT_VECTOR_HPP_INCLUDED
#define BIT_VECTOR_HPP_INCLUDED
#include <cstddef>//std::size_t
#include <cstdint>//std::uint64_t
#include <vector>
std::size_t constexpr bits_per_uint64_t=64;

  //The following bit_vector implementation
  //was copy&pasted from:
  /*
From: Michael Marcin <mike.marcin@gmail.com>
Newsgroups: gmane.comp.lib.boost.devel
Subject: Re: interest in structure of arrays container?
Date: Sun, 30 Oct 2016 18:45:41 -0500
   */
  //and then minor modifications made.
struct bit_vector 
: private std::vector<std::uint64_t>
{
  using super_t=std::vector<std::uint64_t>;
  using super_t::data;
  using super_t::begin;
  using super_t::end;
  
  bit_vector():num_bits_m(0){}  
  void resize( std::size_t num_bits, bool value = false )
  {
      std::size_t num_blocks 
        = (num_bits + (bits_per_uint64_t-1) ) 
        / bits_per_uint64_t;
      constexpr std::uint64_t block_none = 0;
      constexpr std::uint64_t block_all = ~block_none;
      super_t::resize( num_blocks, value ? block_all : block_none );
      num_bits_m = num_bits;
  }
  std::size_t size() const { return num_bits_m; }
  std::size_t vec_size() const { return super_t::size(); }
private:  
  std::size_t num_bits_m;
};   
#endif//BIT_VECTOR_HPP_INCLUDED
