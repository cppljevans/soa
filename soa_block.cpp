//Purpose:
//  Create allocator for single block as described by item:
//    * sizeof...(Ts) allocations could be a single large block 
//  from post:
//    http://lists.boost.org/Archives/boost/2016/10/231136.php
//============================================================
#include <iostream>
#include <array>
#include <cassert>
#include <cmath>
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
  #define TRACE_NXT_OFFSET
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
   *  Ts[vec_size]...
   *  Last offset(at sizeTs) is for the char buffer size.
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
  template
  < typename... Ts
  >
  struct
soa_block
  /**@brief
   *  An *outline* of how to put all arrays in contiguous block
   *  of memory (the my_storage member variable below).
   *  Of course, all the other member functions for
   *  accessing the "internal subarrays" would need to be added.
   */
  {
  private:
    std::size_t my_vec_size;
    using offsets_t=decltype(vec_offsets<Ts...>(1));
    offsets_t my_offsets;
    char*my_storage;
  public:
    soa_block(std::size_t a_vec_size)
    : my_vec_size(a_vec_size)
    , my_offsets(vec_offsets<Ts...>(a_vec_size))
    , my_storage(new char[my_offsets.back()])
      //IIRC, the alignment of the return from new
      //is large enough to accommodate any alignment.
      //Hence, my_offsets.front() can be 0, which
      //is what is returned by vec_offsets.
    {  
      //***TODO***
      //initialize elements of Ts[a_vec_size]...
      //using some sort of placement new code here.
    }
    void resize(std::size_t a_vec_size)
    { offsets_t resized_offsets(vec_offsets<Ts...>(a_vec_size));
      char* resized_storage=new char[resized_offsets.back()];
      //***TODO***
      //copy my_storage to resized_storage,
      //Then, initialize resized_storage that has *not*
      //been copied to, if needed, using, as with
      //CTOR, some sort of placement new code.
      my_offsets=resized_offsets;
      delete[] my_storage;
      my_storage=resized_storage;
    }
  };
  
  template
  < unsigned Id
  >
  struct 
  alignas
  ( 1<<Id%7 
    //==pow(2,Id%7).
    //The c++ standard, IIRC, says all alignments are power of 2.
  ) 
udt //some User Defined Type.
  {
    char c[2*(Id+1)];//Take up some space.
  };
  template
  < typename... Ts
  >
  void
test_offsets(std::size_t vec_size)
  {
    std::cout<<":vec_size="<<vec_size<<":sizeof...(Ts)="<<sizeof...(Ts)<<"\n";
    auto offsets=vec_offsets<Ts...>(vec_size);
    for( auto offset: offsets)
    {
      std::cout<<"offset="<<offset<<"\n";
    }
  }  
  int 
main()
  {
    std::size_t vec_size=5;
    test_offsets<udt<0>,udt<0>>(vec_size);
    test_offsets<udt<0>,udt<3>>(vec_size);
    test_offsets<udt<3>,udt<0>>(vec_size);
    soa_block<udt<3>,udt<0>> soa_v(vec_size);
    soa_v.resize(vec_size*2);
    return 0;
  }  
