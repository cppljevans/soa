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
    std::size_t vec_size=5;
    soa_block<udt<3>,udt<0>> soa_v(vec_size);
    auto n=soa_v.vec_size();
    unsigned i=0;
    for(auto* p=soa_v.begin<0>(); i<n; ++p,++i)
      std::cout<<"p["<<i<<"]="<<*p<<"\n";
    i=0;
    for(auto* p=soa_v.begin<1>(); i<n; ++p,++i)
      std::cout<<"p["<<i<<"]="<<*p<<"\n";
    soa_v.resize(vec_size*2);
    return 0;
  }  
