#ifndef ENUM_SEQUENCE_HPP_INCLUDED
#define ENUM_SEQUENCE_HPP_INCLUDED
//Purpose:
//   create enum_sequence that's like 
//   std::integer_sequence except for enumerations.
#include <utility>//std::integer_sequence
#include <type_traits>//std::underlying_type
  template
  < typename Enumeration
  , Enumeration... Enumerators
  >
  struct
enum_sequence
  {};  
namespace detail_enum_sequence
{ //helper struct's for use by make_enum_sequence. 
    template
    < typename Enumeration
    , typename Indices
    >
    struct
  cast_enum_seq
    ;
    template
    < typename Enumeration
    , typename EnumUnderly
    , EnumUnderly... Indices
    >
    struct
  cast_enum_seq
    < Enumeration
    , std::integer_sequence
      < EnumUnderly
      , Indices...
      >
    >
    {
        using
      type 
        = enum_sequence<Enumeration, Enumeration(Indices)...>
        ;
    };
    template
    < typename Enumeration
    , typename EnumUnderly
    , Enumeration Last
    >
    struct
  cast_enum
    : cast_enum_seq
      < Enumeration
      , std::make_integer_sequence
        < EnumUnderly
        , EnumUnderly(Last)
        >
      >
    {
    };
}//exit detail_enum_sequence namespace  
  template
  < typename Enumeration
  , Enumeration Last
  >
  using
make_enum_sequence=typename
  detail_enum_sequence::cast_enum
  < Enumeration
  , std::underlying_type_t<Enumeration>
  , Last
  >::type
  ;
#endif//ENUM_SEQUENCE_HPP_INCLUDED
