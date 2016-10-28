#ifndef SOA_BLOCK_HPP_INCLUDED
#define SOA_BLOCK_HPP_INCLUDED
#include <cstddef>//std::size_t
#include <boost/align/aligned_alloc.hpp>//aligned_alloc, aligned_free
#include <tuple>
  template
  < std::size_t Index
  , typename T
  >
  struct
soa_vec
  {
  private:
        static  
      T* 
    cast(char*storage)
      {  return static_cast<T*>(static_cast<void*>(storage));
      }
        static
      T const* 
    cast(char const* storage)
      {  return static_cast<T const*>(static_cast<void const*>(storage));
      }
  public:
      T*
    ptr_at(char* my_storage, std::size_t index)
      {
        return cast(my_storage)+index;
      }
      T const*
    ptr_at(char const* my_storage, std::size_t index)const
      {
        return cast(my_storage)+index;
      }
      void
    construct_default
      ( char* my_storage
      , std::size_t my_size=1
      )
      {
        T* begin=cast(my_storage);
        T* end=begin+my_size;
        for(T* now=begin; now!=end; ++now)
          new(now) T;
      }
      void
    construct_copy
      ( char* my_storage
      , char const* other_storage
      , std::size_t size=1
      )
      {
        T* to_begin=cast(my_storage);
        T* to_end=to_begin+size;
        T const*fr_now=cast(other_storage);
        for
          ( T* to_now=to_begin
          ; to_now!=to_end
          ; ++to_now,++fr_now
          )
          new(to_now) T(*fr_now);
      }
      void
    destruct(char* my_storage, std::size_t size=1)
      {
        T* begin=cast(my_storage);
        T* end=begin+size;
        for(T* now=begin; now!=end; ++now)
          now->~T();
      }
      void
    print( std::ostream&sout, char const* my_storage, std::size_t size=1)const
      {
        T const* begin=cast(my_storage);
        T const* end=begin+size;
        unsigned i=0;
        for(T const* now=begin; now!=end; ++now,++i)
          sout<<"  "<<"soa_vec["<<i<<"]="<<*now<<"\n";
      }
  };
  template
  < typename Indices
  , typename... Ts
  >
  struct
soa_impl
  ;
#include "vec_offsets.hpp"
#include <utility>//std::integer_sequence
  template
  < std::size_t... Indices
  , typename... Ts
  >
  struct
soa_impl
  < std::index_sequence< Indices...>
  , Ts...
  >
  : soa_vec<Indices,typename get_type<Ts>::type>... 
  {
  private:
    std::size_t my_vec_size;
    using offsets_t=decltype(vec_offsets<Ts...>(1));
    offsets_t my_offsets;
    void*my_storage;
      template
      < std::size_t Index
      , typename T
      >
        static
      soa_vec<Index,T>*
    get_self
      ( soa_vec<Index,T>*self
      )
      { return self;
      }
      template
      < std::size_t Index
      , typename T
      >
        static
      soa_vec<Index,T> const*
    get_self
      ( soa_vec<Index,T>const*self
      )
      { return self;
      }
      template
      < std::size_t Index
      >
      auto*
    get_vec()
      { return get_self<Index>(this);
      }
      template
      < std::size_t Index
      >
      auto const*
    get_vec()const
      { return get_self<Index>(this);
      }
      char*
    storage_char()
      { return static_cast<char*>(my_storage);
      }
      char const*
    storage_char()const
      { return static_cast<char const*>(my_storage);
      }
      char*
    storage_at(std::size_t index)
      { return storage_char()+my_offsets[index];
      }
      char const*
    storage_at(std::size_t index)const
      { return storage_char()+my_offsets[index];
      }
      template
      < std::size_t Index
      >
      void
    construct_default
      (
      )
      {
        get_vec<Index>()->construct_default(storage_at(Index),my_vec_size);
      }
      template
      < std::size_t Index
      >
      void
    construct_copy
      ( soa_impl const& other_impl
      )
      {
        get_vec<Index>()->construct_copy
        ( storage_at(Index)
        , other_impl.storage_at(Index)
        , my_vec_size
        );
      }
      template
      < std::size_t Index
      >
      void
    destruct()
      {
        get_vec<Index>()->destruct(storage_at(Index),my_vec_size);
      }
      template
      < std::size_t Index
      >
        static
      auto
    mk_index()
      {
        return std::integral_constant<std::size_t,Index>{};
      }
      template
      < typename Fun
      >
        static
      void
    for_each_index(Fun& fun)
     /**@brief
      *   Do fun(mk_index<Indices>())...
      */
      {
        using swallow = int[]; // guaranties left to right order
        (void) //ignore result of following CTOR call.
        swallow
        { 0 //provide something, in case sizeof...(Indices)==0.
        , ( fun(mk_index<Indices>())
          , 0 //ignore result of above fun call and use this value.
          )...
        };
      }
      void*
    my_alloc()const
      { 
        std::size_t const l_align=my_offsets.align_all();
        std::size_t const l_size=my_offsets.size_all();
        void*result=
          boost::alignment::aligned_alloc
          ( l_align
          , l_size
          );
      #if 0
        std::cout<<__func__
          <<":l_align="<<l_align
          <<":l_size="<<l_size
          <<":result="<<result
          <<"\n";
      #endif
        return result;
      }
      void
    my_free()
      { 
        boost::alignment::aligned_free(my_storage);
      }
  protected:
      template
      < std::size_t Index
      >
      void
    print_index(std::ostream& sout)const
      {
        sout<<"soa_block<"<< Index <<">=\n";
        get_vec<Index>()->print(sout,storage_at(Index),my_vec_size);
      }
    soa_impl
      ( std::size_t a_vec_size
      )
      : my_vec_size
        ( a_vec_size
        )
      , my_offsets
        ( vec_offsets<Ts...>
          ( a_vec_size
          )
        )
      , my_storage
        ( my_alloc()
        )
      { 
        auto fun=[this](auto index)
          { //std::cout<<"index()="<<index()<<"\n";
            this->template construct_default<index()>();
          };
        for_each_index(fun);
      }
    soa_impl(soa_impl const& a_other)
      : my_vec_size(a_other.my_vec_size)
      , my_offsets(a_other.my_offsets)
      , my_storage
        ( my_alloc()
        )
      { 
        auto fun = [this,&a_other](auto index)
          { this->template construct_copy<index()>
            ( a_other
            );
          };
        for_each_index(fun);
      }
    ~soa_impl()
      {
        if(my_storage)
        {
          auto fun=[this](auto index)
            { //std::cout<<__func__<<":index()="<<index()<<"\n";
              this->template destruct<index()>();
            };
          for_each_index(fun);
          my_free();
        }
      }
      void
    swap(soa_impl& other)
      {
        std::swap(my_offsets , other.my_offsets);
        std::swap(my_storage , other.my_storage);
        std::swap(my_vec_size, other.my_vec_size);
      }
  public:
      std::size_t
    vec_size()const
      { return my_vec_size;
      }
      soa_impl&
    operator=(soa_impl const& other)
      { 
        //***TODO***
        //  This is probably not the fastest method ;(
        soa_impl copy(other);
        auto fun=[this](auto index)
          { //std::cout<<"index()="<<index()<<"\n";
            this->template destruct<index()>();
          };
        for_each_index(fun);
        swap(copy);
        std::free(copy.my_storage);
        copy.my_storage=0;
        return *this;
      }
      void 
    resize(std::size_t a_vec_size)
      { 
        //***TODO***
        //  This is probably not the fastest method ;(
        soa_impl copy(a_vec_size);
        this->operator=(copy);
      }
      template
      < std::size_t Index
      >
      auto*
    begin()
      {
        return get_vec<Index>()->ptr_at(storage_at(Index),0);
      }
      template
      < std::size_t Index
      >
      auto*
    end()
      {
        return get_vec<Index>()->ptr_at(storage_at(Index),my_vec_size);
      }
      template
      < std::size_t Index
      >
      auto const*
    begin()const
      {
        return get_vec<Index>()->ptr_at(storage_at(Index),0);
      }
      template
      < std::size_t Index
      >
      auto const*
    end()const
      {
        return get_vec<Index>()->ptr_at(storage_at(Index),my_vec_size);
      }
      auto
    begin_all()
      {
        return std::make_tuple(this->template begin<Indices>()...);
      }
      auto
    begin_all()const
      {
        return std::make_tuple(this->template begin<Indices>()...);
      }
      void
    print
      ( std::ostream& sout
      )const
      {
        auto fun=[this,&sout](auto index)
          { this->template print_index<index()>(sout);
          };
        for_each_index(fun);
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
  : public soa_impl
    < std::index_sequence_for<Ts...>
    , Ts...
    >
  {
  public:
      using 
    super_t=
      soa_impl
      < std::index_sequence_for<Ts...>
      , Ts...
      >;
    using super_t::vec_size;
    using super_t::resize;
    
    soa_block
      ( std::size_t a_vec_size
      )
      : super_t
        ( a_vec_size
        )
      {  
      }
        friend
      std::ostream&
    operator<<
      ( std::ostream& sout
      , soa_block const& x
      )
      {
        x.print(sout);
        return sout;
      }
  };
#endif//SOA_BLOCK_HPP_INCLUDED
