#include "bit_vector.hpp"
#include <iostream>
#include <bitset>
struct as_bit_vector_val
{
  std::bitset<bits_per_block> val;
  using from_t=bit_vector::value_type;
  
  as_bit_vector_val(from_t from_v):val(from_v)
    {}
      friend
    std::ostream&
  operator<<
    ( std::ostream&sout
    , as_bit_vector_val const& x
    )
    {
      sout<<x.val;
      return sout;
    }
};
std::ostream& operator<<(std::ostream& sout, bit_vector const& bv)
{
  std::size_t bv_size=bv.size();
  std::size_t vec_size=bv.vec_size();
  std::size_t bits=vec_size*bits_per_block;
  sout
    <<":bv_size="<<bv_size
    <<":vec_size="<<vec_size
    <<":bits="<<bits
    <<"\n"
    ;
  for(unsigned i=0; i<vec_size; ++i)
  {
    sout<<"bit_vector["<<i<<"]=\n";
    sout<<as_bit_vector_val(bv[i])<<"\n";
  }
  return sout;
}
#include <cassert>
#include "sse.hpp"
using sse_vec=sse_float[sse_width];
std::ostream& operator<<(std::ostream& sout, __m128 const& m128)
{
  sse_vec v;
  _mm_store_ps(v,m128);
  sout<<"[";
  for(unsigned i=0; i<sse_width; ++i)
    sout<<i<<":"<<v[i]<<" ";
  sout<<"]";
  return sout;
}
std::ostream& operator<<(std::ostream& sout, sse_vec const& v)
{
  for(unsigned i=0; i<sse_width; ++i)
    sout<<i<<"="<<v[i]<<" ";
  sout<<"\n";
  return sout;
}

#include <string>
#include <iomanip>
unsigned indentation=0;
struct ind
{
  static constexpr unsigned margin=2;
  ind(){ indentation+=margin;}
  ~ind(){ indentation-=margin;}
      friend 
    std::ostream&
  operator<<
    ( std::ostream&sout
    , ind const&
    )
    {
      sout<<std::setw(indentation)<<"";
      return sout;
    }
};  

int main()
{
  std::cout.imbue(std::locale(""));//for thousands separator.
  std::cout<<"sizeof(float)="<<sizeof(float)<<"\n";
  std::cout<<"sizeof(__m128)="<<sizeof(__m128)<<"\n"; 
  std::cout<<"sizeof(uint64_t)="<<sizeof(uint64_t)<<"\n"; 
  std::cout<<"sizeof(std::bitset<bits_per_block>)="<<sizeof(std::bitset<bits_per_block>)<<"\n";
  bit_vector bv;
  std::size_t n=bits_per_block;
    auto 
  update= [](unsigned resize=0)
    { //The following code essentially copy&pasted from
      //the benchmark code for the "SSE_opt" method.
      //IOW, the code starting at line 565 of
      //  http://codepad.org/IbdVcdq8
      //on 2016-10-31.
      //The purpose is to only check to assure block_ptr
      //does not point outside of alive.data.
      std::cout<<"resize="<<resize<<"\n";
      bit_vector alive;
      alive.resize(resize);
      std::cout<<"alive start=\n"<<alive<<"\n";
      std::size_t const n_array //was n in orig. code.
        =alive.size();
      if(n_array%bits_per_block != 0) {
        std::cout<<"***alive.size()="<<n_array<<" must be multiple of "<<bits_per_block<<"\n";
        std::cout<<"***  update aborted.\n";
        return;
      }
      __m128 zero = _mm_setzero_ps();
      std::cout<<"__m128 zero="<<zero<<"\n";
        struct
      e_data
        {
          e_data()
            : val_
              { sse_float(1)
              , sse_float(1)
              , sse_float(1)
              , sse_float(1)
              }
              //Set all values to positive so that in the
              //prints below, it's easy to see how
              //the bitvalues change from the orignal 0's
              //(indicating value is not 0)
              //to 1's 
              //(indicating value is 0).
            , offset_(0)
              //only used to illustrate if offset
              //from 0 is aligned to sse_width.
            {}
            auto
          operator+(std::size_t increment)
            { 
              offset_+=increment;
              return val_;
            }
            bool
          is_aligned()const
            {
              return offset_%sse_width == 0;
            }
        private:
          sse_vec val_;
          std::size_t offset_;
        };
        alignas(sse_align) e_data 
      e_ptr//=energy.data() in original
        ;
        auto 
      e_m128 =
        _mm_load_ps
        ( e_ptr+0
        );
      std::cout<<"__m128 e_m128="<<e_m128<<"\n";
    //{{emulate orig. code
        auto* 
      block_ptr = alive.data()
        //=alive.data() in orig. code
        ;
        using 
      block_t = std::remove_reference<decltype(*block_ptr)>::type
        //= uint64_t in orig. code
        ;
      auto* block_end = alive.end();//absent in orig. code
      for 
      ( size_t i_array //was i in orig. code
        = 0
      ; ( i_array < n_array 
          //requires n_array%sse_width==0 to work correctly.
          //if not, the i_array>=n_array in
          //inner loop.
       && block_ptr<block_end
          //absent in orig. code
        )
      ; ++block_ptr
        //Absent in orig. code.
        //Instead, there was *block_ptr++ = block_v
        //just after the, in orig. code, do...while statement.
      ) 
      {
        ind a_ind;
        std::cout<<a_ind<<"i_array="<<i_array<<"\n";
        block_t block_v = 0;
        for //was a do...while loop in orig. code
        ( size_t i_bit=0//absent in orig. code
        ; 
          ( i_bit < bits_per_block
          && i_array < n_array//unneeded if bits_per_block%sse_width==0.
          )//was while ( i % bits_per_block != 0 ) in orig. code
        ; ( i_array+=sse_width
          , i_bit+=sse_width//absent in orig. code
          )
        )
        {
          ind a_ind;
          std::cout<<a_ind<<"i_bit="<<i_bit<<"\n";
          auto e_i = e_ptr + i_array;
          std::cout<<a_ind<<"e_i="<<e_i<<"\n";
          assert(e_ptr.is_aligned());//absent in orig. code
          auto load_e=
            _mm_load_ps
            ( e_i
            );
          std::cout<<a_ind<<"load_e="<<load_e<<"\n";
          auto is_positive = 
            _mm_cmpgt_ps
            ( load_e
            , zero
            );
          int mask = _mm_movemask_ps(is_positive);
          std::cout<<a_ind<<"mask="<<mask<<"\n";
            block_t 
          shiftee
            ( mask
              //Was _mm_movemask_ps (.. e_i) in orig code,
              //which sets on/off most signficant sse_width bits in orig code.
              //References:
              //  https://msdn.microsoft.com/en-us/library/4490ys29(v=vs.90).aspx
            )
            ;
          std::cout<<a_ind<<"shiftee="<<as_bit_vector_val(shiftee)<<"\n";
          block_t bitor_right = shiftee << i_bit;
          std::cout<<a_ind
            <<"bitor_right["<<i_array<<"]="
            <<as_bit_vector_val(bitor_right)
            <<"\n";
          block_v |= bitor_right;
          std::cout<<a_ind
            <<"block["<<i_array<<"]="
            <<as_bit_vector_val(block_v)
            <<"\n";
        }
        assert(block_ptr<block_end);
        *block_ptr=block_v;
        std::cout<<a_ind<<"*block_ptr="<<as_bit_vector_val(*block_ptr)<<"\n";
      }
      std::cout<<"alive finish=\n"<<alive<<"\n";
    //}}emulate orig. code
    };//update
  update(0);
  update(1);
  update(n);
  update(2*n);
  update(n+1);
  return 0;
}
