#ifndef VEC_OFFSETS_HPP_INCLUDED
#define VEC_OFFSETS_HPP_INCLUDED
#include <array>//and std::make_tuple, std::tuple 
#include <boost/align/is_aligned.hpp>
#include <boost/align/align_up.hpp>
namespace vec_offsets_alignment=boost::alignment;
  template
  < typename T
  , std::size_t Align=alignof(T) 
    //Requested alignment for T.
    //Must be (some power of 2)*alignof(T)
  >
  struct
type_align
  /**@brief
   *  Allow overalignment, Align for T
   *  to be specified in template args 
   *  to vec_offsets.
   */
  {
    using type=T;
    static constexpr std::size_t align=Align;
  };
  template
  < typename T
  >
  struct
get_align
  {
    static constexpr std::size_t align=alignof(T);
  };
  template
  < typename T
  >
  struct
get_type
  {
    using type=T;
  };
  template
  < typename T
  , std::size_t Align
  >
  struct
get_align
  < type_align<T,Align>
  >
  {
    static constexpr std::size_t align=Align;
  };
  template
  < typename T
  , std::size_t Align
  >
  struct
get_type
  < type_align<T,Align>
  >
  {
    using type=T;
  };
  template
  < std::size_t N
  >
  using
over_aligns_t= 
  std::array
  < std::size_t
  , N
  >
  ;
  template
  < std::size_t N
  >
  struct
vec_offsets_t
  : std::array
    < std::size_t
    , N+2
    >
  {
    vec_offsets_t()
      {}
      static constexpr std::size_t
    i_size=N;
      static constexpr std::size_t
    i_align=i_size+1;
    std::size_t size_all()const
      { return this->operator[](i_size);}
    std::size_t align_all()const
      { return this->operator[](i_align);}
  };
  template
  < typename... TypeAlign 
  >
    constexpr
  auto
vec_offsets
  ( std::size_t vec_size
  )
  /**@brief
   *  Calculate offsets in a char buffer for storing 
   *  suitably aligned get_type<TypeAlign>::type... elements.
   *  Next to last value is is for the char buffer size.
   *  Last value is maximum of get_align<TypeAlign>::align...
   */
  {
    std::size_t const num_types=sizeof...(TypeAlign);
    using offsets_t=vec_offsets_t<num_types>;
    offsets_t offsets_v;
    std::size_t sizes[num_types]={(vec_size*sizeof(typename get_type<TypeAlign>::type))...};
    std::size_t aligns[num_types]={get_align<TypeAlign>::align...};
    std::size_t i_type=0;
    std::size_t max_align=aligns[i_type];
    offsets_v[i_type]=0;
      auto 
    nxt_offset=[&](std::size_t align_i)
      {
        offsets_v[i_type]=vec_offsets_alignment::align_up
          ( offsets_v[i_type-1]+sizes[i_type-1]
          , align_i
          );
      };
    for
      ( ++i_type
      ; i_type<num_types
      ; ++i_type
      )
    {
      nxt_offset(aligns[i_type]);
      max_align=std::max(max_align,aligns[i_type]);
    }
    nxt_offset(max_align);
    offsets_v[offsets_t::i_align]=max_align;
    return offsets_v;
  }
#endif//VEC_OFFSETS_HPP_INCLUDED
