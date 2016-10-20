#ifndef UDT_HPP_INCLUDED
#define UDT_HPP_INCLUDED
#include <algorithm>//std::max
#include <iostream>
  using
udt_id_t=unsigned
  ;  
  udt_id_t
udt_count=0
 /**@brief
  *  Count of udt's created in order to provide
  *  a unique my_id for each one created.
  */
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
    , std::alignment_of<udt_id_t>::value
      //account for my_id alignment.
    )
  ) 
udt //some User Defined Type.
  {
  private:
    udt_id_t my_id;
    char c[2*(Id+1)];//Take up some space.
  public:
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
      udt&
    operator=(udt const&)
      {
        //do nothing to preserve uniqueness of my_id.
        return *this;
      }
    ~udt()
      {
        std::cout<<"DTOR:"<<*this<<"\n";
      }
  };
#endif//UDT_HPP_INCLUDED
