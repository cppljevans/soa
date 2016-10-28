#ifndef UDT_HPP_INCLUDED
#define UDT_HPP_INCLUDED
#include <algorithm>//std::max
#include <set>
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
  std::multiset<udt_id_t>
is_live
  ;
  std::ostream&
operator<<
  ( std::ostream& sout
  , std::multiset<udt_id_t> const& x
  )
  {
    for(auto e: x)
    {
      sout<<e<<"\n";
    }
    return sout;
  }  
  struct
udt_id_base
  {
  private:
    udt_id_t my_id;
  public:
      udt_id_t
    id()const
      { return my_id;}
    udt_id_base()
      : my_id(udt_count++)
      {
        is_live.insert(my_id);
      }
    udt_id_base(udt_id_base&)
      : my_id(udt_count++)
      {
        is_live.insert(my_id);
      }
    ~udt_id_base()
      {
        is_live.erase(is_live.find(my_id));
      }
      udt_id_base&
    operator=(udt_id_base const&)
      { //don't copy my_id to preserve uniqueness.
        return *this;
      }
  };
  template
  < unsigned TypeId//type identifier
  >
  struct 
  alignas
  ( std::max
    ( std::size_t(1<<TypeId%7)
      //==pow(2,Id%7).  
      //The 7 is somewhat arbitrary.
      //It's purpose is just to limit size of alignment.
      //The c++ standard, IIRC, says all alignments are power of 2.
    , std::alignment_of<udt_id_t>::value
      //account for my_id alignment.
    )
  ) 
udt //some User Defined Type.
  : private udt_id_base
  {
  private:
    udt_id_t my_val;//used to test if copy done.
    char c[2*(TypeId+1)];//Take up some space.
  public:
        friend
      std::ostream&
    operator<<
      ( std::ostream& sout
      , udt const& x
      )
      {
        return 
          sout
          <<"udt<"<<TypeId<<">"
          <<":my_id="<<x.id()
          <<":my_val="<<x.my_val
          <<":sizeof="<<sizeof(udt)
          <<":alignof="<<alignof(udt)
          <<":this="<<(void const*)&x
          ;
      }
    udt()
      : my_val(this->udt_id_base::id())
      {
        std::cout<<"CTOR(default):*this="<<*this<<"\n";
      }
    udt(udt const& x)
      : my_val(x.my_val)
      {
        std::cout<<"CTOR(copy):"<<*this<<"\n";
      }
      udt&
    operator=(udt const& x)
      {
        my_val=x.val;//to show copy actually done.
        return *this;
      }
    ~udt()
      {
        std::cout<<"DTOR:"<<*this<<"\n";
      }
  };
#endif//UDT_HPP_INCLUDED
