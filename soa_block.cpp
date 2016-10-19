//Purpose:
//  Create allocator for single block as described by item:
//    * sizeof...(Ts) allocations could be a single large block 
//  from post:
//    http://lists.boost.org/Archives/boost/2016/10/231136.php
//============================================================
#include <iostream>
#include <array>
#include <cassert>
#include <type_traits>
#include <algorithm>
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
  < std::size_t Index
  , typename T
  >
  struct
soa_vec
  {
  public:
      void
    construct(char* storage, std::size_t size=1)
      {
        T* begin=cast(storage);
        T* end=begin+size;
        for(T* now=begin; now!=end; ++now)
          new(now) T;
      }
      void
    destruct(char* storage, std::size_t size=1)
      {
        T* begin=cast(storage);
        T* end=begin+size;
        for(T* now=begin; now!=end; ++now)
          now->~T();
      }
      T*
    ptr_at(char* storage, std::size_t index)
      {
        return cast(storage)+index;
      }
      T const*
    ptr_at(char const* storage, std::size_t index)const
      {
        return cast(storage)+index;
      }
  private:      
      T* 
    cast(char* storage)
      {  return static_cast<T*>(static_cast<void*>(storage));
      }
      T const* 
    cast(char const* storage)
      {  return static_cast<T const*>(static_cast<void const*>(storage));
      }
  };
  template
  < typename Indices
  , typename... Ts
  >
  struct
soa_impl
  ;
  template
  < std::size_t... Indices
  , typename... Ts
  >
  struct
soa_impl
  < std::integer_sequence< std::size_t, Indices...>
  , Ts...
  >
  : soa_vec<Indices,Ts>... 
  {
  private:
    std::size_t my_vec_size;
    using offsets_t=decltype(vec_offsets<Ts...>(1));
    offsets_t my_offsets;
    char*my_storage;
      template
      < std::size_t Index
      , typename T
      >
        static
      soa_vec<Index,T>&
    get_self
      ( soa_vec<Index,T>&self
      )
      { return self;
      }
      template
      < std::size_t Index
      , typename T
      >
        static
      soa_vec<Index,T> const&
    get_self
      ( soa_vec<Index,T>const&self
      )
      { return self;
      }
      template
      < std::size_t Index
      >
      auto&
    get_vec()
      { return get_self<Index>(*this);
      }
      template
      < std::size_t Index
      >
      auto const&
    get_vec()const
      { return get_self<Index>(*this);
      }
      char*
    storage_at(std::size_t index)
      { return my_storage+my_offsets[index];
      }
      char const*
    storage_at(std::size_t index)const
      { return my_storage+my_offsets[index];
      }
      template
      < std::size_t Index
      >
      void
    construct
      (
      )
      {
        get_vec<Index>().construct(storage_at(Index),my_vec_size);
      }
      template
      < std::size_t Index
      >
      void
    destruct
      (
      )
      {
        get_vec<Index>().destruct(storage_at(Index),my_vec_size);
      }
  protected:
    soa_impl(std::size_t a_vec_size)
      : my_vec_size(a_vec_size)
      , my_offsets(vec_offsets<Ts...>(a_vec_size))
      , my_storage(new char[my_offsets.back()])
        //IIRC, the alignment of the return from new
        //is large enough to accommodate any alignment.
        //Hence, my_offsets.front() can be 0, which
        //is what is returned by vec_offsets.
      { 
        using swallow = int[]; // guaranties left to right order
        (void)swallow{0, (void(this->template construct<Indices>()), 0)...};
      }
      void 
    resize(std::size_t a_vec_size)
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
    ~soa_impl()
      {
        using swallow = int[]; // guaranties left to right order
        (void)swallow{0, (void(this->template destruct<Indices>()), 0)...};
        delete[] my_storage;
      }
  public:
      template
      < std::size_t Index
      >
      auto*
    begin()
      {
        return get_vec<Index>().ptr_at(storage_at(Index),0);
      }
      template
      < std::size_t Index
      >
      auto*
    end()
      {
        return get_vec<Index>().ptr_at(storage_at(Index),my_vec_size);
      }
      template
      < std::size_t Index
      >
      auto const*
    begin()const
      {
        return get_vec<Index>().ptr_at(storage_at(Index),0);
      }
      template
      < std::size_t Index
      >
      auto const*
    end()const
      {
        return get_vec<Index>().ptr_at(storage_at(Index),my_vec_size);
      }
  };
  template
  < typename... Ts
  >
  struct
soa_block
  /**@brief
   *  Implementation of the single block of storage for
   *  structure of arrays mentioned in //Purpose comment
   *  at top of this file.
   */
  : public soa_impl< std::index_sequence_for<Ts...>, Ts...>
  {
  public:
    using super_t=soa_impl< std::index_sequence_for<Ts...>, Ts...>;
    soa_block(std::size_t a_vec_size)
    : super_t(a_vec_size)
    {  
    }
    void resize(std::size_t a_vec_size)
    { this->super_t::resize(a_vec_size);
    }
  };
  
  unsigned
udt_count=0
  ;
  template
  < unsigned Id
  >
  struct 
  alignas
  ( std::max
    ( std::size_t(1<<Id%7)
      //==pow(2,Id%7).
      //The c++ standard, IIRC, says all alignments are power of 2.
    , std::alignment_of<unsigned>::value
      //account for my_id alignment.
    )
  ) 
udt //some User Defined Type.
  {
    unsigned my_id;
    char c[2*(Id+1)];//Take up some space.
        friend
      std::ostream&
    operator<<
      ( std::ostream& sout
      , udt const& x
      )
      {
        return sout<<"udt<"<<Id<<">:my_id="<<x.my_id;
      }
    udt()
      : my_id(udt_count++)
      {
        std::cout<<"CTOR(default):"<<*this<<"\n";
      }
    udt(udt const&)
      : my_id(udt_count++)
      {
        std::cout<<"CTOR(copy):"<<*this<<"\n";
      }
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
    unsigned i=0;
    for(auto* p=soa_v.begin<0>(); p< soa_v.end<0>(); ++p,++i)
      std::cout<<"p["<<i<<"]="<<*p<<"\n";
    i=0;
    for(auto* p=soa_v.begin<1>(); p< soa_v.end<1>(); ++p,++i)
      std::cout<<"p["<<i<<"]="<<*p<<"\n";
    soa_v.resize(vec_size*2);
    return 0;
  }  
