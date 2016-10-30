#ifndef INDEX_TYPE_BUF_HPP_INCLUDED
#define INDEX_TYPE_BUF_HPP_INCLUDED
#include <cstddef>//std::size_t
#include <iostream>
  template
  < std::size_t Index
  , typename T
  >
  struct
index_type_buf
  /**@brief
   *  Provide methods for a T[N] stored in buffer.
   *  N is defined by subclass.
   */
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
          sout<<"  "<<"index_type_buf["<<i<<"]="<<*now<<"\n";
      }
  };
#endif//INDEX_TYPE_BUF_HPP_INCLUDED
