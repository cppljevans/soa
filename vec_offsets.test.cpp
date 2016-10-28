#include "vec_offsets.hpp"
#include <iostream>
#include <string>
#include <cassert>
  template
  < std::size_t OverAlign//==some power of 2.
  , typename... Ts
  >
  void
test_offsets
  ( std::size_t vec_size
  )
  {
    constexpr std::size_t num_types=sizeof...(Ts);
    std::cout<<":vec_size="<<vec_size<<":sizeof...(Ts)="<<num_types<<"\n";
    auto offsets=
      vec_offsets
      < type_align
        < Ts
        , alignof(Ts)*OverAlign
        >...
      >
      ( vec_size
      );
    std::size_t aligns[]={alignof(Ts)...};
    unsigned n=offsets.size()-2;
    std::cout<<"offsets.size()-2="<<n<<"\n";
    assert(n==num_types);
    unsigned i=0;
    for( i=0; i<n; ++i)
    {
      std::size_t offset=offsets[i];
      std::cout<<"offset["<<i<<"]="<<offset<<"\n";
      std::cout<<"alignment["<<i<<"]="<<aligns[i]<<"\n";
      std::size_t over_alignment=aligns[i]*OverAlign;
      std::cout<<"over_alignment["<<i<<"]="<<over_alignment<<"\n";
      bool aligned=vec_offsets_alignment::is_aligned(offsets[i],aligns[i]);
      std::cout<<"aligned["<<i<<"]="<<aligned<<"\n";
      bool over_aligned=vec_offsets_alignment::is_aligned(offsets[i],over_alignment);
      std::cout<<"over_aligned["<<i<<"]="<<over_aligned<<"\n";
      assert(aligned);
      assert(over_aligned);
    }
    std::size_t buf_size=offsets[n];
    std::cout<<"buf_size="<<buf_size<<"\n";
    std::size_t max_align=std::max(alignof(Ts)...);
    std::cout<<"max_alignment="<<max_align<<"\n";
    bool max_aligned=vec_offsets_alignment::is_aligned(buf_size,max_align);
    std::cout<<"max_aligned="<<max_aligned<<"\n";
    if(!max_aligned)
    {
      unsigned unaligned=0;
      std::string indent="***";
      for( i=0; i<n; ++i)
      {
        std::size_t offset=buf_size+offsets[i];
        std::cout<<indent<<"buf_size+offset["<<i<<"]="<<offset<<"\n";
        std::cout<<indent<<"alignment["<<i<<"]="<<aligns[i]<<"\n";
        bool aligned=vec_offsets_alignment::is_aligned(offset,aligns[i]);
        if(!aligned) ++unaligned;
        std::cout<<indent<<"aligned["<<i<<"]="<<aligned<<"\n";
      }
      std::cout<<indent<<"unaligned="<<unaligned<<"\n";
      assert(max_aligned);
    }
  }  
#include "udt.hpp"
#include <iomanip>
  template
  < std::size_t OverAlignLog2
  >
  void
test_overaligned
  ()
  {
    std::cout<<"===OverAlignLog2="<<OverAlignLog2<<"\n";
    constexpr std::size_t over_align=std::size_t(1)<<OverAlignLog2;//==some power of 2.
    std::cout<<"===over_align="<<over_align<<"\n";
    std::cout<<std::boolalpha;
    std::size_t vec_size=3;
    test_offsets<over_align,udt<0>,udt<0>>(vec_size);
    test_offsets<over_align,udt<0>,udt<3>>(vec_size);
    test_offsets<over_align,udt<3>,udt<0>>(vec_size);
    test_offsets<over_align,udt<4>,udt<2>>(vec_size);
  }
  int 
main()
  {
    test_overaligned<0>();
    test_overaligned<1>();
    test_overaligned<2>();
    return 0;
  }
