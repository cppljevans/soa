//Purpose:
//  Create allocator for single block as described by item:
//    * sizeof...(Ts) allocations could be a single large block 
//  from post:
//    http://lists.boost.org/Archives/boost/2016/10/231136.php
//============================================================
#include <iostream>
#include "soa_block.hpp"
#include "udt.hpp"  
  int 
main()
{
  {
    std::size_t vec_size=5;
    using soa_t=soa_block<udt<3>,udt<0>>;
    soa_t soa_v0(vec_size);
    std::cout<<"soa_v0=\n"<<soa_v0<<"\n";
    soa_t soa_v1(soa_v0);
    std::cout<<"soa_v1=\n"<<soa_v1<<"\n";
    soa_t soa_v2(vec_size+2);
    std::cout<<"soa_v2=\n"<<soa_v2<<"\n";
    soa_v1=soa_v2;
    std::cout<<"soa_v1=soa_v2\n"<<soa_v1<<"\n";
    soa_v0.resize(vec_size+3);
    std::cout<<"soa_v0.resize=\n"<<soa_v0<<"\n";
  }
  if(is_live.size() !=0)
  {
    std::cout<<"***ERROR***: is_live should be empty!\n";
    std::cout<<"is_live=\n"<<is_live<<".\n";
  }
  return 0;
}  
