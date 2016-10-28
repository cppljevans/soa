//Purpose:
//  Compare methods of Structure Of Arrays data structure.
//OriginalSource:
//  On 2016-10-19 by cppljevans
//    WHAT: downloaded original source from:
//      http://codepad.org/eol6auRN
//RESULT:
/*
/tmp/build/clangxx3_8_pkg/clang/struct_of_arrays/work/soa_compare.benchmark.exe
particle_count=1,000,000
minimum duration=5.18095

comparitive performance table:

method   rel_duration
________ ______________
Block    1
StdArray 1.00159
Flat     1.00874
SoA      1.0325
AoS      1.39591

Compilation finished at Mon Oct 24 08:26:48
 */
//=============================
#define NDEBUG //disable assert's.
#if defined(__GNUC__)
  #include <emmintrin.h>
#else
  #include <mmintrin.h>
#endif//for __m128 type and _mm_* functions.
#include <vector>
#include <string>
#include <chrono>
#include <iostream>
#include <locale>
#include <random>
#include <algorithm>
#include <array>
//#define HAVE_GOON_BIT_VECTOR
#ifdef HAVE_GOON_BIT_VECTOR
  #include <goon/bit_vector.hpp>
#endif
#include <boost/align/aligned_allocator.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>

using namespace std;
using boost::alignment::aligned_allocator;
#ifdef HAVE_GOON_BIT_VECTOR
  using goon::bit_vector;
#else
  using bit_vector = std::vector<char>;
#endif

struct float2_t {
    float x, y;
};

struct float3_t {
    float x, y, z;

    friend float3_t operator-( const float3_t & lhs, const float3_t & rhs ) {
        return float3_t{ lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
    }
    friend float3_t operator+( const float3_t & lhs, const float3_t & rhs ) {
        return float3_t{ lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
    }
    friend float3_t operator*( const float3_t & lhs, float rhs ) {
        return float3_t{ lhs.x * rhs, lhs.y * rhs, lhs.z * rhs };
    }

    float3_t& operator+=( const float3_t & rhs ) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
};

struct float4_t {
    float x, y, z, w;
};

struct particle_nest_t {
    float3_t position;
    float3_t velocity;
    float3_t acceleration;
    float2_t size;
    float4_t color;
    float energy;
    bool alive;
};

struct particle_nest_enum
{
  enum particle_e {
      position,
      velocity,
      acceleration,
      size,
      color,
      energy,
      alive
  };
};

struct particle_flat_enum
{
  enum particle_e {
      position_x,
      position_y,
      position_z,
      velocity_x,
      velocity_y,
      velocity_z,
      acceleration_x,
      acceleration_y,
      acceleration_z,
      size,
      color,
      energy,
      alive
  };
};
uniform_real_distribution<float> pdist( -10.f, 10.f );
uniform_real_distribution<float> vxzdist( -5.f, 5.f );
uniform_real_distribution<float> vydist( 4.f, 20.f );
uniform_real_distribution<float> adist( -1.f, 1.f );
uniform_real_distribution<float> sdist( 0.1f, 5.f );
uniform_real_distribution<float> cdist(  0.f, 1.f );
uniform_real_distribution<float> edist( 10.f, 1000.f );

constexpr float gravity = -9.8f;
constexpr float dt = 1.0f;
//#define USE_SMALL_PARTICLE_COUNT
#ifdef USE_SMALL_PARTICLE_COUNT
constexpr size_t particle_count = 1024;
constexpr size_t frames = 1000;
#else
constexpr size_t particle_count = 15625 * 64;
constexpr size_t frames = 1000;
#endif

enum
method_enum
  { AoS
  , SoA_vec
  , SoA_flat
  , SoA_stdArr
  , SoA_block
  , SSE_block
  , SSE_vec
  , SSEopt_vec
  , LFA
  , method_last
  };

  template
  < method_enum Method
  >
struct emitter_t
  ;
  template<>
struct emitter_t<AoS> {
    static constexpr char const*name(){return "AoS";}
    static constexpr bool defined=true;

    vector<particle_nest_t> particles;

    void generate( size_t n, mt19937 & rng ) {
        particles.resize( n );


        for ( size_t i = 0; i < n; ++i ) {
            particles[i].position = float3_t{ pdist( rng ), pdist( rng ), pdist( rng ) };
            particles[i].velocity = float3_t{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            particles[i].acceleration = float3_t{ adist( rng ), adist( rng ), adist( rng ) };
            particles[i].size = float2_t{ sdist( rng ), sdist( rng ) };
            particles[i].color = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            particles[i].energy = edist( rng );
            particles[i].alive = true;
        }
    }

    void update() {
        for ( auto & prt : particles ) {
            auto & p = prt.position;
            auto & v = prt.velocity;
            auto & a = prt.acceleration;
            auto & e = prt.energy;
            v += (a * dt);
            v.y += gravity * dt;

            p += (v * dt);
            e -= dt;
            if ( e <= 0 ) {
                prt.alive = false;
            }
        }
    }
};

  template<>
struct emitter_t<SoA_vec> {
    static constexpr char const*name(){return "SoA_vec";}
    static constexpr bool defined=true;

    vector<float3_t> position;
    vector<float3_t> velocity;
    vector<float3_t> acceleration;
    vector<float2_t> size;
    vector<float4_t> color;
    vector<float>  energy;
    vector<char>   alive;

    void generate( size_t n, mt19937 & rng ) {
        position.resize( n );
        velocity.resize( n );
        acceleration.resize( n );
        size.resize( n );
        color.resize( n );
        energy.resize( n );
        alive.resize( n, true );

        for ( size_t i = 0; i < n; ++i ) {
            position[i] = float3_t{ pdist( rng ), pdist( rng ), pdist( rng ) };
            velocity[i] = float3_t{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            acceleration[i] = float3_t{ adist( rng ), adist( rng ), adist( rng ) };
            size[i] = float2_t{ sdist( rng ), sdist( rng ) };
            color[i] = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            energy[i] = edist( rng );
        }
    }

    void update() {

        size_t n = position.size();
        auto ps = position.data();
        auto vs = velocity.data();
        auto as = acceleration.data();
        auto es = energy.data();
        auto als = alive.data();
        for ( size_t i = 0; i < n; ++i ) {
            auto & p = ps[i];
            auto & v = vs[i];
            auto & a = as[i];
            auto & e = es[i];
            v += (a * dt);
            v.y += gravity * dt;

            p += (v * dt);
            e -= dt;
            if ( e <= 0 ) {
                als[i] = false;
            }
        }
    }

};

  template<>
struct emitter_t<SoA_flat> {
    static constexpr char const*name(){return "SoA_flat";}
    static constexpr bool defined=true;

    char * data;
    size_t capacity;

    void free() {
        delete[] data;
    }

    ~emitter_t() {
        free();
    }

    void resize( size_t n ) {
        capacity = n;
        data = new char[
            sizeof( float3_t )*n +
            sizeof( float3_t )*n +
            sizeof( float3_t )*n +
            sizeof( float2_t )*n +
            sizeof( float4_t )*n +
            sizeof( float )*n +
            sizeof( char )*n
        ];
    }

    float3_t* get_position() 
      { return reinterpret_cast<float3_t*>(data); }
    float3_t* get_velocity() 
      { return reinterpret_cast<float3_t*>(data + capacity * offsetof(particle_nest_t, velocity)); }
    float3_t* get_acceleration()
      { return reinterpret_cast<float3_t*>(data + capacity * offsetof(particle_nest_t, acceleration)); }
    float2_t* get_size()
      { return reinterpret_cast<float2_t*>(data + capacity * offsetof(particle_nest_t, size)); }
    float4_t* get_color()
      { return reinterpret_cast<float4_t*>(data + capacity * offsetof(particle_nest_t, color)); }
    float* get_energy()
      { return reinterpret_cast<float*>(data + capacity * offsetof(particle_nest_t, energy)); }
    char* get_alive()
      { return reinterpret_cast<char*>(data + capacity * offsetof(particle_nest_t, alive)); }

    void generate( size_t n, mt19937 & rng ) {
        resize( n );

        for ( size_t i = 0; i < n; ++i ) {
            get_position()[i] = float3_t{ pdist( rng ), pdist( rng ), pdist( rng ) };
            get_velocity()[i] = float3_t{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            get_acceleration()[i] = float3_t{ adist( rng ), adist( rng ), adist( rng ) };
            get_size()[i] = float2_t{ sdist( rng ), sdist( rng ) };
            get_color()[i] = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            get_energy()[i] = edist( rng );
            get_alive()[i] = true;
        }
    }

    void update() {

        size_t n = capacity;
        auto ps = get_position();
        auto vs = get_velocity();
        auto as = get_acceleration();
        auto es = get_energy();
        auto als = get_alive();
        for ( size_t i = 0; i < n; ++i ) {
            auto & p = ps[i];
            auto & v = vs[i];
            auto & a = as[i];
            auto & e = es[i];
            v += (a * dt);
            v.y += gravity * dt;

            p += (v * dt);
            e -= dt;
            if ( e <= 0 ) {
                als[i] = false;
            }
        }
    }
};

template< typename T >
using soa_array = array<T, particle_count>;
  template<>
struct emitter_t<SoA_stdArr> 
: public particle_nest_enum
{
    static constexpr char const*name(){return "SoA_stdArr";}
    static constexpr bool defined=true;

    typedef tuple<
        soa_array<float3_t>,
        soa_array<float3_t>,
        soa_array<float3_t>,
        soa_array<float2_t>,
        soa_array<float4_t>,
        soa_array<float>,
        soa_array<char> > particles_t;

    unique_ptr<particles_t> particles;

      template<particle_e Index>
    auto& get(){ return std::get<Index>(*particles);}
    
    void resize( size_t n) {
        particles.reset( new particles_t );
    }
    
    std::size_t particles_size() {
        return get<(particle_e)0>().size();
    }
    
    void generate( size_t n, mt19937 & rng ) {
        resize( n );
        
        for ( size_t i = 0; i < n; ++i ) {
            get<position>()[i] = float3_t{ pdist( rng ), pdist( rng ), pdist( rng ) };
            get<velocity>()[i] = float3_t{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            get<acceleration>()[i] = float3_t{ adist( rng ), adist( rng ), adist( rng ) };
            get<size>()[i] = float2_t{ sdist( rng ), sdist( rng ) };
            get<color>()[i] = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            get<energy>()[i] = edist( rng );
            get<alive>()[i] = true;
        }
    }

    void update() {
        size_t n = particles_size();

        for ( size_t i = 0; i < n; ++i ) {
            auto & p = get<position>()[i];
            auto & v = get<velocity>()[i];
            auto & a = get<acceleration>()[i];
            auto & e = get<energy>()[i];
            v += (a * dt);
            v.y += gravity * dt;

            p += (v * dt);
            e -= dt;
            if ( e <= 0 ) {
                get<alive>()[i] = false;
            }
        }
    }
};

#include "soa_block.hpp"
  template<>
struct emitter_t<SoA_block>
/**@brief
 *  Pretty much cut&past from above
 *    soa_emitter_static_t
 *  but use soa_block instead of soa_array.
 */
: public particle_nest_enum
{
    static constexpr char const*name(){return "SoA_block";}
    static constexpr bool defined=true;
    
    typedef soa_block<
        float3_t,
        float3_t,
        float3_t,
        float2_t,
        float4_t,
        float,
        bool> data_t;

    data_t data;

    emitter_t()
      : data(particle_count)
      {}

    void generate( size_t n, mt19937 & rng ) {
        auto begins_v=data.begin_all();
        for ( size_t i = 0; i < n; ++i ) {
            get<position>(begins_v)[i] = float3_t{ pdist( rng ), pdist( rng ), pdist( rng ) };
            get<velocity>(begins_v)[i] = float3_t{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            get<acceleration>(begins_v)[i] = float3_t{ adist( rng ), adist( rng ), adist( rng ) };
            get<size>(begins_v)[i] = float2_t{ sdist( rng ), sdist( rng ) };
            get<color>(begins_v)[i] = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            get<energy>(begins_v)[i] = edist( rng );
            get<alive>(begins_v)[i] = true;
        }
    }

    void update() {
        auto begins_v=data.begin_all();
        size_t n = particle_count;
        for ( size_t i = 0; i < n; ++i ) {
            auto & p = get<position>(begins_v)[i];
            auto & v = get<velocity>(begins_v)[i];
            auto & a = get<acceleration>(begins_v)[i];
            auto & e = get<energy>(begins_v)[i];
            v += (a * dt);
            v.y += gravity * dt;

            p += (v * dt);
            e -= dt;
            if ( e <= 0 ) {
                get<alive>(begins_v)[i] = false;
            }
        }
    }
};

#include <libflatarray/flat_array.hpp>
class particle_lfa_t
{
public:
    float position[3];
    float velocity[3];
    float acceleration[3];
    float size[2];
    float color[4];
    float energy;
    float alive;
};

LIBFLATARRAY_REGISTER_SOA(
    particle_lfa_t,
    ((float)(position)(3))
    ((float)(velocity)(3))
    ((float)(acceleration)(3))
    ((float)(size)(2))
    ((float)(color)(4))
    ((float)(energy))
    ((float)(alive))
)

  template<>
struct emitter_t<LFA> 
{
    static constexpr char const*name(){return "LFA";}
    static constexpr bool defined=true;
    LibFlatArray::soa_vector<particle_lfa_t> particles;

    void generate( size_t n, mt19937 & rng ) 
    {
        particles.resize( n );

        // The SoA layout is fixed at compile time. Multiple layouts
        // are available via templates. callback() dispatches at
        // runtime to the correct template instantiation:
        particles.callback([n, &rng](auto particle){
                for (; particle.index() < n; ++particle ) {
                    particle.position()[0] = pdist( rng );
                    particle.position()[1] = pdist( rng );
                    particle.position()[2] = pdist( rng );
                    particle.velocity()[0] = vxzdist( rng );
                    particle.velocity()[1] = vxzdist( rng );
                    particle.velocity()[2] = vxzdist( rng );
                    particle.acceleration()[0] = adist( rng );
                    particle.acceleration()[1] = adist( rng );
                    particle.acceleration()[2] = adist( rng );
                    particle.size()[0] = sdist( rng );
                    particle.size()[1] = sdist( rng );
                    particle.color()[0] = cdist( rng );
                    particle.color()[1] = cdist( rng );
                    particle.color()[2] = cdist( rng );
                    particle.color()[3] = cdist( rng );
                    particle.energy() = edist( rng );
                    particle.alive() = 1;
                }
            });
    }

    void update() 
    {
      using LibFlatArray::any;
      using LibFlatArray::get;
  
      long n = particles.size();
      long hits = 0;
      particles.callback
      ( [n, &hits](auto particle)
        {
          // short_vec behaves much like an ordinary float, but
          // maps all computation to vector intrinsics. The ISA
          // (SSE, AVX, AVX512...) is chosen automatically at
          // compile time.
          typedef LibFlatArray::short_vec<float, 16> Float;
    
          // The loop peeler handles left-over loop iterations
          // when the number of particles is not a multiple of
          // the vector arity...
          LibFlatArray::loop_peeler<Float>
          ( &particle.index()
          , n
          , [&particle, n, &hits](auto new_float, long *i, long end) 
            {
              // ...by switching arities:
              typedef decltype(new_float) Float;
              Float dt2 = dt;
              for (; particle.index() < end; particle += Float::ARITY) 
              {
                // vector loads are issued by passing a pointer to the c-tor:
                Float v0 = Float(&particle.velocity()[0]) + (Float(&particle.acceleration()[0]) * dt2);
                Float v1 = Float(&particle.velocity()[1]) + (Float(&particle.acceleration()[1]) * dt2);
                Float v2 = Float(&particle.velocity()[2]) + (Float(&particle.acceleration()[2]) * dt2);
                v1 += gravity * dt2;
  
                &particle.velocity()[0] << v0;
                &particle.velocity()[1] << v1;
                &particle.velocity()[2] << v2;
  
                Float p0 = Float(&particle.position()[0]) + (Float(&particle.velocity()[0]) * dt2);
                Float p1 = Float(&particle.position()[1]) + (Float(&particle.velocity()[1]) * dt2);
                Float p2 = Float(&particle.position()[2]) + (Float(&particle.velocity()[2]) * dt2);
  
                &particle.position()[0] << p0;
                &particle.position()[1] << p1;
                &particle.position()[2] << p2;
  
                Float e = Float(&particle.energy()) - dt2;
                &particle.energy() << e;
  
                // initialization from scalar value
                // broadcasts the value to all vector
                // lanes:
                Float alive = 1;
                // Comparison creates a bit-mask which can
                // be used to selectively set values in
                // the target register(s):
                alive.blend((e <= 0.0f), 0.0);
                &particle.alive() << alive;
              }
            }
          );
        }
      );
    }
};

std::size_t constexpr sse_align=16;
  template<>
struct emitter_t<SSE_block>
: public particle_flat_enum
{
    static constexpr char const*name(){return "SSE_block";}
    static constexpr bool defined=true;
      soa_block
      <
        type_align<float,sse_align>,// position_x;
        type_align<float,sse_align>,// position_y;
        type_align<float,sse_align>,// position_z;
        type_align<float,sse_align>,// velocity_x;
        type_align<float,sse_align>,// velocity_y;
        type_align<float,sse_align>,// velocity_z;
        type_align<float,sse_align>,// acceleration_x;
        type_align<float,sse_align>,// acceleration_y;
        type_align<float,sse_align>,// acceleration_z;
        float2_t,// size;
        float4_t,// color;
        type_align<float,sse_align>,// energy;
        char// alive;
      > particles;
      
      template<particle_e Index>
      auto* get(){ return particles.begin<Index>();}
      template<particle_e Index>
      auto* data(){ return particles.begin<Index>();}
      
    void resize( size_t n) {
        particles.resize(n);
    }
    
    std::size_t particles_size() {
        return particles.vec_size();
    }
    
    void generate( size_t n, mt19937 & rng ) {
        resize( n );

        for ( size_t i = 0; i < n; ++i ) {
            get<position_x>()[i] = pdist( rng );
            get<position_y>()[i] = pdist( rng );
            get<position_z>()[i] = pdist( rng );
            get<velocity_x>()[i] = vxzdist( rng );
            get<velocity_y>()[i] = vxzdist( rng );
            get<velocity_z>()[i] = vxzdist( rng );
            get<acceleration_x>()[i] = adist( rng );
            get<acceleration_y>()[i] = adist( rng );
            get<acceleration_z>()[i] = adist( rng );
            get<size>()[i] = float2_t{ sdist( rng ), sdist( rng ) };
            get<color>()[i] = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            get<energy>()[i] = edist( rng );
        }
    }
    
    void update() {
        size_t n = particles_size();
        __m128 vx, vy, vz;
        __m128 t = _mm_set1_ps( dt );
        __m128 g = _mm_set1_ps( gravity * dt );
        __m128 zero = _mm_setzero_ps();
        for ( size_t i = 0; i < n; i += 4 ) {
            vx = _mm_add_ps( _mm_load_ps( data<velocity_x>()+i ), _mm_mul_ps( t, _mm_load_ps( data<acceleration_x>()+i )));
            vy = _mm_add_ps( _mm_load_ps( data<velocity_y>()+i ), _mm_mul_ps( t, _mm_load_ps( data<acceleration_y>()+i )));
            vz = _mm_add_ps( _mm_load_ps( data<velocity_z>()+i ), _mm_mul_ps( t, _mm_load_ps( data<acceleration_z>()+i )));

            vy = _mm_add_ps( vy, g );

            _mm_store_ps( data<position_x>()+i, _mm_add_ps(_mm_load_ps(data<position_x>()+i), _mm_mul_ps(vx, t)));
            _mm_store_ps( data<position_y>()+i, _mm_add_ps(_mm_load_ps(data<position_y>()+i), _mm_mul_ps(vy, t)));
            _mm_store_ps( data<position_z>()+i, _mm_add_ps(_mm_load_ps(data<position_z>()+i), _mm_mul_ps(vz, t)));

            _mm_store_ps( data<velocity_x>()+i, vx );
            _mm_store_ps( data<velocity_y>()+i, vy );
            _mm_store_ps( data<velocity_z>()+i, vz );

            _mm_store_ps( data<energy>()+i, _mm_sub_ps(_mm_load_ps(data<energy>()+i), t));

            auto a = _mm_movemask_ps( _mm_cmple_ps( _mm_load_ps( data<energy>()+i ), zero ));
            for ( int j = 0; j < 4; ++j ) {
                get<alive>()[i+j] = (a & (1 << j));
            }
        }
    }
};

template< typename T >
using sse_vector = vector<T, aligned_allocator<T,sse_align> >;

  template<>
struct emitter_t<SSE_vec> 
: public particle_flat_enum
{
    static constexpr char const*name(){return "SSE_vec";}
    static constexpr bool defined=true;
      tuple
      <
        sse_vector<float>,// position_x;
        sse_vector<float>,// position_y;
        sse_vector<float>,// position_z;
        sse_vector<float>,// velocity_x;
        sse_vector<float>,// velocity_y;
        sse_vector<float>,// velocity_z;
        sse_vector<float>,// acceleration_x;
        sse_vector<float>,// acceleration_y;
        sse_vector<float>,// acceleration_z;
        vector<float2_t>,// size;
        vector<float4_t>,// color;
        sse_vector<float>,// energy;
        vector<char>// alive;
      > particles;
      
      template<particle_e Index>
      auto& get(){ return std::get<Index>(particles);}
      template<particle_e Index>
      auto* data(){ return std::get<Index>(particles).data();}
      
      void resize( size_t n) {
        get<position_x>().resize( n );
        get<position_y>().resize( n );
        get<position_z>().resize( n );
        get<velocity_x>().resize( n );
        get<velocity_y>().resize( n );
        get<velocity_z>().resize( n );
        get<acceleration_x>().resize( n );
        get<acceleration_y>().resize( n );
        get<acceleration_z>().resize( n );
        get<size>().resize( n );
        get<color>().resize( n );
        get<energy>().resize( n );
        get<alive>().resize( n, true );
    }
    
    void generate( size_t n, mt19937 & rng ) {
        resize( n );

        for ( size_t i = 0; i < n; ++i ) {
            get<position_x>()[i] = pdist( rng );
            get<position_y>()[i] = pdist( rng );
            get<position_z>()[i] = pdist( rng );
            get<velocity_x>()[i] = vxzdist( rng );
            get<velocity_y>()[i] = vxzdist( rng );
            get<velocity_z>()[i] = vxzdist( rng );
            get<acceleration_x>()[i] = adist( rng );
            get<acceleration_y>()[i] = adist( rng );
            get<acceleration_z>()[i] = adist( rng );
            get<size>()[i] = float2_t{ sdist( rng ), sdist( rng ) };
            get<color>()[i] = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            get<energy>()[i] = edist( rng );
        }
    }

    std::size_t particles_size() {
        return get<position_x>().size();
    }
    
    void update() {
        size_t n = particles_size();
        __m128 vx, vy, vz;
        __m128 t = _mm_set1_ps( dt );
        __m128 g = _mm_set1_ps( gravity * dt );
        __m128 zero = _mm_setzero_ps();
        for ( size_t i = 0; i < n; i += 4 ) {
            vx = _mm_add_ps( _mm_load_ps( data<velocity_x>()+i ), _mm_mul_ps( t, _mm_load_ps( data<acceleration_x>()+i )));
            vy = _mm_add_ps( _mm_load_ps( data<velocity_y>()+i ), _mm_mul_ps( t, _mm_load_ps( data<acceleration_y>()+i )));
            vz = _mm_add_ps( _mm_load_ps( data<velocity_z>()+i ), _mm_mul_ps( t, _mm_load_ps( data<acceleration_z>()+i )));

            vy = _mm_add_ps( vy, g );

            _mm_store_ps( data<position_x>()+i, _mm_add_ps(_mm_load_ps(data<position_x>()+i), _mm_mul_ps(vx, t)));
            _mm_store_ps( data<position_y>()+i, _mm_add_ps(_mm_load_ps(data<position_y>()+i), _mm_mul_ps(vy, t)));
            _mm_store_ps( data<position_z>()+i, _mm_add_ps(_mm_load_ps(data<position_z>()+i), _mm_mul_ps(vz, t)));

            _mm_store_ps( data<velocity_x>()+i, vx );
            _mm_store_ps( data<velocity_y>()+i, vy );
            _mm_store_ps( data<velocity_z>()+i, vz );

            _mm_store_ps( data<energy>()+i, _mm_sub_ps(_mm_load_ps(data<energy>()+i), t));

            auto a = _mm_movemask_ps( _mm_cmple_ps( _mm_load_ps( data<energy>()+i ), zero ));
            for ( int j = 0; j < 4; ++j ) {
                get<alive>()[i+j] = (a & (1 << j));
            }
        }
    }
};

  template<>
struct emitter_t<SSEopt_vec> {
    static constexpr char const*name(){return "SSEopt_vec";}
    static constexpr bool defined=true;

    sse_vector<float> position_x;
    sse_vector<float> position_y;
    sse_vector<float> position_z;
    sse_vector<float> velocity_x;
    sse_vector<float> velocity_y;
    sse_vector<float> velocity_z;
    sse_vector<float> acceleration_x;
    sse_vector<float> acceleration_y;
    sse_vector<float> acceleration_z;
    vector<float2_t> size;
    vector<float4_t> color;
    sse_vector<float>  energy;
    bit_vector alive;

    void generate( size_t n, mt19937 & rng ) {
        position_x.resize( n );
        position_y.resize( n );
        position_z.resize( n );
        velocity_x.resize( n );
        velocity_y.resize( n );
        velocity_z.resize( n );
        acceleration_x.resize( n );
        acceleration_y.resize( n );
        acceleration_z.resize( n );
        size.resize( n );
        color.resize( n );
        energy.resize( n );
        alive.resize( n, true );

        for ( size_t i = 0; i < n; ++i ) {
            position_x[i] = pdist( rng );
            position_y[i] = pdist( rng );
            position_z[i] = pdist( rng );
            velocity_x[i] = vxzdist( rng );
            velocity_y[i] = vxzdist( rng );
            velocity_z[i] = vxzdist( rng );
            acceleration_x[i] = adist( rng );
            acceleration_y[i] = adist( rng );
            acceleration_z[i] = adist( rng );
            size[i] = float2_t{ sdist( rng ), sdist( rng ) };
            color[i] = float4_t{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            energy[i] = edist( rng );
        }
    }

    void update() {
        size_t n = position_x.size();
        __m128 vx, vy, vz;
        __m128 t = _mm_set1_ps( dt );
        __m128 g = _mm_set1_ps( gravity * dt );
        for ( size_t i = 0; i < n; i += 4 ) {
            vx = _mm_add_ps( _mm_load_ps( velocity_x.data()+i ), _mm_mul_ps( t, _mm_load_ps( acceleration_x.data()+i )));
            vy = _mm_add_ps( _mm_load_ps( velocity_y.data()+i ), _mm_mul_ps( t, _mm_load_ps( acceleration_y.data()+i )));
            vz = _mm_add_ps( _mm_load_ps( velocity_z.data()+i ), _mm_mul_ps( t, _mm_load_ps( acceleration_z.data()+i )));

            vy = _mm_add_ps( vy, g );

            _mm_store_ps( position_x.data()+i, _mm_add_ps( _mm_load_ps( position_x.data()+i ), _mm_mul_ps( vx, t )));
            _mm_store_ps( position_y.data()+i, _mm_add_ps( _mm_load_ps( position_y.data()+i ), _mm_mul_ps( vy, t )));
            _mm_store_ps( position_z.data()+i, _mm_add_ps( _mm_load_ps( position_z.data()+i ), _mm_mul_ps( vz, t )));

            _mm_store_ps( velocity_x.data()+i, vx );
            _mm_store_ps( velocity_y.data()+i, vy );
            _mm_store_ps( velocity_z.data()+i, vz );
        }
        __m128 zero = _mm_setzero_ps();
        uint64_t *block_ptr = (uint64_t*)alive.data();
        auto e_ptr = energy.data();
        for ( size_t i = 0; i < n; ) {
            uint64_t block = 0;
            do {
                _mm_store_ps( e_ptr + i, _mm_sub_ps( _mm_load_ps( e_ptr + i ), t ));
                block |= uint64_t( _mm_movemask_ps( _mm_cmple_ps( _mm_load_ps( e_ptr + i ), zero ))) << (i % 64);
                i += 4;
            } while ( i % 64 != 0 );
            *block_ptr++ = block;
        }
    }
};

#include <utility>
  using
dur_t=double;
  struct
run_result_t
  {
    const char* method;//name of method;
    dur_t dur_v;//time taken by method;
  };
#include <limits>
  constexpr auto 
dur_max=std::numeric_limits<dur_t>::max();
  template< typename emitter_t >
  run_result_t
run_test()
  {
    using clock_t = chrono::high_resolution_clock;

    const auto seed_val = mt19937::default_seed;
    mt19937 rng( seed_val );

    dur_t diff=dur_max;
    bool defined=emitter_t::defined;
    if(defined)
    {
      emitter_t emitter;
      emitter.generate( particle_count, rng );
  
      auto start = clock_t::now();
      
      for ( size_t i = 0; i < frames; ++i ) 
      {
          emitter.update();
      }
      auto finish = clock_t::now();
      chrono::duration<dur_t> dur(finish-start);
      diff=dur.count();
    }
    return run_result_t{emitter_t::name(),diff};
  }
#include <iomanip>
  template
  < std::size_t N
  >
  void
print_results
  ( std::ostream&sout
  , std::vector<run_result_t>& run_result_v
  )
  {
     auto compare=[](run_result_t const& i, run_result_t const& j)
       { return (i.dur_v < j.dur_v);
       };
     std::sort(run_result_v.begin(),run_result_v.end(),compare);
     dur_t dur_min=run_result_v[0].dur_v;
     std::cout<<"minimum duration="<<dur_min<<"\n\n";
     std::cout<<"comparitive performance table:\n\n";
     unsigned const ncol=3;
     std::string const header[ncol]={"method   ","duration","rel_dur"};
     unsigned long wcol[ncol];
     sout<<std::left;
     for(unsigned i=0; i<ncol; ++i)
     {
       wcol[i]=header[i].size()+2;
       sout
         <<std::setw(wcol[i])<<header[i]
         <<" ";
     }
     sout<<"\n";
     sout<<std::setfill('_');
     for(unsigned i=0; i<ncol; ++i)
     {
       sout
         <<std::setw(wcol[i])<<""
         <<" ";
     }
     sout<<"\n";
     sout<<std::setfill(' ');
     int const prec=6;
     for(unsigned i=0; i<N; ++i)
     {
       sout
         <<std::setw(wcol[0])<<run_result_v[i].method
         <<" ";
       auto dur_i=run_result_v[i].dur_v;
       bool is_defined = (dur_i < dur_max*.9);
       if(is_defined)
       {
         sout
           <<std::setw(wcol[1])<<dur_i
           <<" "
           <<std::setw(wcol[2])<<std::setprecision(prec)<<(dur_i/dur_min)
           ;
       }
       else
       {
         sout
           <<std::setw(wcol[2])<<"undefined"
           <<" "
           <<std::setw(wcol[2])<<"undefined"
           ;
       }
       sout
         <<"\n";
     }
  }
#include "enum_sequence.hpp"
  template
  < method_enum... Methods
  >
  void
run_tests
  ( enum_sequence<method_enum,Methods...>
  )
  {  
     std::vector<run_result_t> run_result_v={run_test<emitter_t<Methods>>()...};
     print_results<sizeof...(Methods)>(std::cout,run_result_v);
  }
int main()
  {
    cout.imbue(std::locale(""));//for thousands separator.
    cout << "particle_count="<< particle_count << std::endl;
    cout << "frames="<< frames << std::endl;
    run_tests
      ( make_enum_sequence
        < method_enum
        , method_last
        >{}
      );
    return 0;
  }
