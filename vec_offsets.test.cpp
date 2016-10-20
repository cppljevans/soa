#include "vec_offsets.hpp"
#include <iostream>
#include <string>
  template
  < typename... Ts
  >
  void
test_offsets(std::size_t vec_size)
  {
    std::cout<<":vec_size="<<vec_size<<":sizeof...(Ts)="<<sizeof...(Ts)<<"\n";
    auto offsets=vec_offsets<Ts...>(vec_size);
    std::size_t aligns[]={alignof(Ts)...};
    unsigned n=offsets.size()-1;
    std::cout<<"sizeof...(Ts)="<<n<<"\n";
    unsigned i=0;
    for( i=0; i<n; ++i)
    {
      std::size_t offset=offsets[i];
      std::cout<<"offset["<<i<<"]="<<offset<<"\n";
      std::cout<<"alignment["<<i<<"]="<<aligns[i]<<"\n";
      bool aligned=offsets[i]%aligns[i]==0;
      std::cout<<"aligned["<<i<<"]="<<aligned<<"\n";
    }
    std::size_t buf_size=offsets[n];
    std::cout<<"buf_size="<<buf_size<<"\n";
    std::size_t max_align=std::max(alignof(Ts)...);
    std::cout<<"max_alignment="<<max_align<<"\n";
    bool max_aligned=buf_size%max_align==0;
    std::cout<<"max_aligned="<<max_aligned<<"\n";
    if(!max_aligned)
    {
      std::string indent="  ";
      for( i=0; i<n; ++i)
      {
        std::size_t offset=buf_size+offsets[i];
        std::cout<<indent<<"buf_size+offset["<<i<<"]="<<offset<<"\n";
        std::cout<<indent<<"alignment["<<i<<"]="<<aligns[i]<<"\n";
        bool aligned=offset%aligns[i]==0;
        std::cout<<indent<<"aligned["<<i<<"]="<<aligned<<"\n";
      }
    }
  }  
#include "udt.hpp"
  int 
main()
  {
    std::size_t vec_size=3;
    test_offsets<udt<0>,udt<0>>(vec_size);
    test_offsets<udt<0>,udt<3>>(vec_size);
    test_offsets<udt<3>,udt<0>>(vec_size);
    test_offsets<udt<4>,udt<2>>(vec_size);
    return 0;
  }
