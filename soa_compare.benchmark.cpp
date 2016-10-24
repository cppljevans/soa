//ChangeLog:
//  On 2016-10-19 by cppljevans
//    WHAT: downloaded original source from:
//      http://codepad.org/eol6auRN
//  On 2016-10-20.1811 by cppljevans:
//    WHAT:
//      added HAVE_* macros
//      and #ifdef HAVE_* ..#endif's
//    WHY:
//      to avoid compile errors due to not having some symbols
//      on my OS.
//  On 2016-10-20.1834 by cppljevans:
//    WHAT:
//      added USE_SMALL_PARTICLE_COUNT macro.
//    WHY:
//      to avoid long run times during development.
//  On 2016-10-20.2137 by cppljevans:
//    WHAT:
//      added struct soa_emitter_block_t
//    WHY:
//      to test performance of soa_block method.
//    RESULT:
//      took about twice as long as soa_emitter_static_t
/*
particle_count=1,000
AoS in 0.0359925 seconds
SoA in 0.0324623 seconds
SoA flat in 0.0325361 seconds
SoA Static in 0.0426253 seconds
SoA block in 0.0805211 seconds
 */
//  On 2016-10-21.0243 by cppljevans:
//    WHAT:
//      Used newly added soa_block.begin_all,
//      and -O2 optimization flag (instead of previous -O0 flag).
//    RESULT:
//      Now, soa_block method fastest(by small amount):
/*
/home/evansl/dwnlds/llvm/3.8/prebuilt/clang+llvm-3.8.0-x86_64-linux-gnu-ubuntu-14.04/bin/clang++ -c -O2 -stdlib=libc++  -std=c++14 -ftemplate-backtrace-limit=0 -fdiagnostics-show-template-tree -fno-elide-type -fmacro-backtrace-limit=0       -DTYPE_AT_IMPL=0   -ftemplate-depth=200  codepad.eol6auRN.cpp -MMD -o /tmp/build/gcc4_9_0/clang/struct_of_arrays/work/codepad.eol6auRN.o
/home/evansl/dwnlds/llvm/3.8/prebuilt/clang+llvm-3.8.0-x86_64-linux-gnu-ubuntu-14.04/bin/clang++ -stdlib=libc++    /tmp/build/gcc4_9_0/clang/struct_of_arrays/work/codepad.eol6auRN.o   -o /tmp/build/gcc4_9_0/clang/struct_of_arrays/work/codepad.eol6auRN.exe
/tmp/build/gcc4_9_0/clang/struct_of_arrays/work/codepad.eol6auRN.exe
particle_count=1,000,000
AoS in 7.18701 seconds
SoA in 5.38178 seconds
SoA flat in 5.19252 seconds
SoA Static in 5.14276 seconds
SoA block in 5.12777 seconds

Compilation finished at Fri Oct 21 02:41:38
 */
//=============================
#include <mmintrin.h>
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

using namespace std;
using boost::alignment::aligned_allocator;
#ifdef HAVE_GOON_BIT_VECTOR
  using goon::bit_vector;
#else
  using bit_vector = std::vector<char>;
#endif

#ifdef HAVE_LIBFLATARRAY
#include <libflatarray/flat_array.hpp>
#endif

struct float2 {
    float x, y;
};

struct float3 {
    float x, y, z;

    friend float3 operator-( const float3 & lhs, const float3 & rhs ) {
        return float3{ lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
    }
    friend float3 operator+( const float3 & lhs, const float3 & rhs ) {
        return float3{ lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
    }
    friend float3 operator*( const float3 & lhs, float rhs ) {
        return float3{ lhs.x * rhs, lhs.y * rhs, lhs.z * rhs };
    }

    float3& operator+=( const float3 & rhs ) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
};

struct float4 {
    float x, y, z, w;
};

struct particle_t {
    float3 position;
    float3 velocity;
    float3 acceleration;
    float2 size;
    float4 color;
    float energy;
    bool alive;
};

#ifdef HAVE_LIBFLATARRAY

class plain_particle_t
{
public:
    float position[3];
    float velocity[3];
    float acceleration[3];
    float size[2];
    float color[4];
    float energy;
    bool alive;
};

LIBFLATARRAY_REGISTER_SOA(
    plain_particle_t,
    ((float)(position)(3))
    ((float)(velocity)(3))
    ((float)(acceleration)(3))
    ((float)(size)(2))
    ((float)(color)(4))
    ((float)(energy))
    ((bool)(alive))
)

#endif

uniform_real_distribution<float> pdist( -10.f, 10.f );
uniform_real_distribution<float> vxzdist( -5.f, 5.f );
uniform_real_distribution<float> vydist( 4.f, 20.f );
uniform_real_distribution<float> adist( -1.f, 1.f );
uniform_real_distribution<float> sdist( 0.1f, 5.f );
uniform_real_distribution<float> cdist(  0.f, 1.f );
uniform_real_distribution<float> edist( 10.f, 1000.f );

constexpr float gravity = -9.8f;
constexpr float dt = 1.0f;
constexpr size_t particle_count =
//#define USE_SMALL_PARTICLE_COUNT
#ifdef USE_SMALL_PARTICLE_COUNT
  1000
#else
  15625 * 64
#endif
  ;

struct aos_emitter_t {
    vector<particle_t> particles;

    void generate( size_t n, mt19937 & rng ) {
        particles.resize( n );


        for ( size_t i = 0; i < n; ++i ) {
            particles[i].position = float3{ pdist( rng ), pdist( rng ), pdist( rng ) };
            particles[i].velocity = float3{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            particles[i].acceleration = float3{ adist( rng ), adist( rng ), adist( rng ) };
            particles[i].size = float2{ sdist( rng ), sdist( rng ) };
            particles[i].color = float4{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
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

struct soa_emitter_t {
    vector<float3> position;
    vector<float3> velocity;
    vector<float3> acceleration;
    vector<float2> size;
    vector<float4> color;
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
            position[i] = float3{ pdist( rng ), pdist( rng ), pdist( rng ) };
            velocity[i] = float3{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            acceleration[i] = float3{ adist( rng ), adist( rng ), adist( rng ) };
            size[i] = float2{ sdist( rng ), sdist( rng ) };
            color[i] = float4{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
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

struct soa_emitter_flat_t {

    char * data;
    size_t capacity;

    ~soa_emitter_flat_t() {
        free();
    }

    void alloc( size_t n ) {
        capacity = n;
        data = new char[
            sizeof( float3 )*n +
            sizeof( float3 )*n +
            sizeof( float3 )*n +
            sizeof( float2 )*n +
            sizeof( float4 )*n +
            sizeof( float )*n +
            sizeof( char )*n
        ];
    }

    void free() {
        delete data;
    }

    float3* get_position()      { return reinterpret_cast<float3*>(data); }
    float3* get_velocity()      { return reinterpret_cast<float3*>(data + capacity * offsetof(particle_t, velocity)); }
    float3* get_acceleration()  { return reinterpret_cast<float3*>(data + capacity * offsetof(particle_t, acceleration)); }
    float2* get_size()          { return reinterpret_cast<float2*>(data + capacity * offsetof(particle_t, size)); }
    float4* get_color()         { return reinterpret_cast<float4*>(data + capacity * offsetof(particle_t, color)); }
    float*  get_energy()        { return reinterpret_cast<float*>(data + capacity * offsetof(particle_t, energy)); }
    char*   get_alive()         { return reinterpret_cast<char*>(data + capacity * offsetof(particle_t, alive)); }

    void generate( size_t n, mt19937 & rng ) {
        alloc( n );

        for ( size_t i = 0; i < n; ++i ) {
            get_position()[i] = float3{ pdist( rng ), pdist( rng ), pdist( rng ) };
            get_velocity()[i] = float3{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            get_acceleration()[i] = float3{ adist( rng ), adist( rng ), adist( rng ) };
            get_size()[i] = float2{ sdist( rng ), sdist( rng ) };
            get_color()[i] = float4{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
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

    enum aspect : std::size_t {
        position,
        velocity,
        acceleration,
        size,
        color,
        energy,
        alive
    };

template< typename T >
using soa_array = array<T, particle_count>;

struct soa_emitter_static_t {
    typedef tuple<
        soa_array<float3>,
        soa_array<float3>,
        soa_array<float3>,
        soa_array<float2>,
        soa_array<float4>,
        soa_array<float>,
        soa_array<char> > data_t;

    unique_ptr<data_t> data;

    void generate( size_t n, mt19937 & rng ) {
        data.reset( new data_t );
        for ( size_t i = 0; i < n; ++i ) {
            get<position>(*data)[i] = float3{ pdist( rng ), pdist( rng ), pdist( rng ) };
            get<velocity>(*data)[i] = float3{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            get<acceleration>(*data)[i] = float3{ adist( rng ), adist( rng ), adist( rng ) };
            get<size>(*data)[i] = float2{ sdist( rng ), sdist( rng ) };
            get<color>(*data)[i] = float4{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            get<energy>(*data)[i] = edist( rng );
            get<alive>(*data)[i] = true;
        }
    }

    void update() {
        size_t n = particle_count;

        for ( size_t i = 0; i < n; ++i ) {
            auto & p = get<position>(*data)[i];
            auto & v = get<velocity>(*data)[i];
            auto & a = get<acceleration>(*data)[i];
            auto & e = get<energy>(*data)[i];
            v += (a * dt);
            v.y += gravity * dt;

            p += (v * dt);
            e -= dt;
            if ( e <= 0 ) {
                get<alive>(*data)[i] = false;
            }
        }
    }
};

#include "soa_block.hpp"
struct soa_emitter_block_t
/**@brief
 *  Pretty much cut&past from above
 *    soa_emitter_static_t
 *  but use soa_block instead of soa_array.
 */
{
    typedef soa_block<
        float3,
        float3,
        float3,
        float2,
        float4,
        float,
        char> data_t;

    data_t data;
    typedef decltype(data.begin_all()) begins_t;
    begins_t const begins_v;

    soa_emitter_block_t()
      : data(particle_count)
      , begins_v(data.begin_all())
      {}

    void generate( size_t n, mt19937 & rng ) {
        for ( size_t i = 0; i < n; ++i ) {
            get<position>(begins_v)[i] = float3{ pdist( rng ), pdist( rng ), pdist( rng ) };
            get<velocity>(begins_v)[i] = float3{ vxzdist( rng ), vydist( rng ), vxzdist( rng ) };
            get<acceleration>(begins_v)[i] = float3{ adist( rng ), adist( rng ), adist( rng ) };
            get<size>(begins_v)[i] = float2{ sdist( rng ), sdist( rng ) };
            get<color>(begins_v)[i] = float4{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            get<energy>(begins_v)[i] = edist( rng );
            get<alive>(begins_v)[i] = true;
        }
    }

    void update() {
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

struct libflatarray_soa_emitter_t {
    LibFlatArray::soa_vector<plain_particle_t> particles;

    void generate( size_t n, mt19937 & rng ) {
        particles.resize( n );

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
                    particle.alive() = true;
                }
            });
    }

    void update() {
        long n = particles.size();
        particles.callback([n](auto particle){
                for (; particle.index() < n; ++particle) {
                    float v0 = particle.velocity()[0] + (particle.acceleration()[0] * dt);
                    float v1 = particle.velocity()[1] + (particle.acceleration()[1] * dt);
                    float v2 = particle.velocity()[2] + (particle.acceleration()[2] * dt);
                    v1 += gravity * dt;

                    particle.velocity()[0] = v0;
                    particle.velocity()[1] = v1;
                    particle.velocity()[2] = v2;

                    float p0 = particle.position()[0] + (particle.velocity()[0] * dt);
                    float p1 = particle.position()[1] + (particle.velocity()[1] * dt);
                    float p2 = particle.position()[2] + (particle.velocity()[2] * dt);

                    particle.position()[0] = p0;
                    particle.position()[1] = p1;
                    particle.position()[2] = p2;

                    float e = particle.energy() - dt;
                    particle.energy() = e;

                    if ( e <= 0 ) {
                        particle.alive() = false;
                    }
                }
            });
    }
};

//#define HAVE__M128
#ifdef HAVE__M128
template< typename T >
using sse_vector = vector<T, aligned_allocator<T,16> >;

struct soa_emitter_sse_t {
    sse_vector<float> position_x;
    sse_vector<float> position_y;
    sse_vector<float> position_z;
    sse_vector<float> velocity_x;
    sse_vector<float> velocity_y;
    sse_vector<float> velocity_z;
    sse_vector<float> acceleration_x;
    sse_vector<float> acceleration_y;
    sse_vector<float> acceleration_z;
    vector<float2> size;
    vector<float4> color;
    sse_vector<float>  energy;
    vector<char>   alive;

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
            size[i] = float2{ sdist( rng ), sdist( rng ) };
            color[i] = float4{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            energy[i] = edist( rng );
        }
    }

    void update() {
        size_t n = position_x.size();

        __m128 vx, vy, vz;
        __m128 t = _mm_set1_ps( dt );
        __m128 g = _mm_set1_ps( gravity * dt );
        __m128 zero = _mm_setzero_ps();
        for ( size_t i = 0; i < n; i += 4 ) {
            vx = _mm_add_ps( _mm_load_ps( velocity_x.data()+i ), _mm_mul_ps( t, _mm_load_ps( acceleration_x.data()+i ) ) );
            vy = _mm_add_ps( _mm_load_ps( velocity_y.data()+i ), _mm_mul_ps( t, _mm_load_ps( acceleration_y.data()+i ) ) );
            vz = _mm_add_ps( _mm_load_ps( velocity_z.data()+i ), _mm_mul_ps( t, _mm_load_ps( acceleration_z.data()+i ) ) );

            vy = _mm_add_ps( vy, g );

            _mm_store_ps( position_x.data()+i, _mm_add_ps(_mm_load_ps(position_x.data()+i), _mm_mul_ps(vx, t)) );
            _mm_store_ps( position_y.data()+i, _mm_add_ps(_mm_load_ps(position_y.data()+i), _mm_mul_ps(vy, t)) );
            _mm_store_ps( position_z.data()+i, _mm_add_ps(_mm_load_ps(position_z.data()+i), _mm_mul_ps(vz, t)) );

            _mm_store_ps( velocity_x.data()+i, vx );
            _mm_store_ps( velocity_y.data()+i, vy );
            _mm_store_ps( velocity_z.data()+i, vz );

            _mm_store_ps( energy.data()+i, _mm_sub_ps(_mm_load_ps(energy.data()+i), t) );

            auto a = _mm_movemask_ps( _mm_cmple_ps( _mm_load_ps( energy.data() + i ), zero ) );
            for ( int j = 0; j < 4; ++j ) {
                alive[i+j] = (a & (1 << j));
            }
        }
    }
};

struct soa_emitter_sse_opt_t {
    sse_vector<float> position_x;
    sse_vector<float> position_y;
    sse_vector<float> position_z;
    sse_vector<float> velocity_x;
    sse_vector<float> velocity_y;
    sse_vector<float> velocity_z;
    sse_vector<float> acceleration_x;
    sse_vector<float> acceleration_y;
    sse_vector<float> acceleration_z;
    vector<float2> size;
    vector<float4> color;
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
            size[i] = float2{ sdist( rng ), sdist( rng ) };
            color[i] = float4{ cdist( rng ), cdist( rng ), cdist( rng ), cdist( rng ) };
            energy[i] = edist( rng );
        }
    }

    void update() {
        size_t n = position_x.size();

        __m128 vx, vy, vz;
        __m128 t = _mm_set1_ps( dt );
        __m128 g = _mm_set1_ps( gravity * dt );

        for ( size_t i = 0; i < n; i += 4 ) {
            vx = _mm_add_ps( _mm_load_ps( velocity_x.data() + i ), _mm_mul_ps( t, _mm_load_ps( acceleration_x.data() + i ) ) );
            vy = _mm_add_ps( _mm_load_ps( velocity_y.data() + i ), _mm_mul_ps( t, _mm_load_ps( acceleration_y.data() + i ) ) );
            vz = _mm_add_ps( _mm_load_ps( velocity_z.data() + i ), _mm_mul_ps( t, _mm_load_ps( acceleration_z.data() + i ) ) );

            vy = _mm_add_ps( vy, g );

            _mm_store_ps( position_x.data() + i, _mm_add_ps( _mm_load_ps( position_x.data() + i ), _mm_mul_ps( vx, t ) ) );
            _mm_store_ps( position_y.data() + i, _mm_add_ps( _mm_load_ps( position_y.data() + i ), _mm_mul_ps( vy, t ) ) );
            _mm_store_ps( position_z.data() + i, _mm_add_ps( _mm_load_ps( position_z.data() + i ), _mm_mul_ps( vz, t ) ) );

            _mm_store_ps( velocity_x.data() + i, vx );
            _mm_store_ps( velocity_y.data() + i, vy );
            _mm_store_ps( velocity_z.data() + i, vz );
        }
        __m128 zero = _mm_setzero_ps();
        uint64_t *block_ptr = (uint64_t*)alive.data();
        auto e_ptr = energy.data();
        for ( size_t i = 0; i < n; ) {
            uint64_t block = 0;
            do {
                _mm_store_ps( e_ptr + i, _mm_sub_ps( _mm_load_ps( e_ptr + i ), t ) );
                block |= uint64_t( _mm_movemask_ps( _mm_cmple_ps( _mm_load_ps( e_ptr + i ), zero ) ) ) << (i % 64);
                i += 4;
            } while ( i % 64 != 0 );
            *block_ptr++ = block;
        }
    }
};
#endif //HAVE__M128

template< typename emitter_t >
chrono::duration<double> run_test() {
    using clock_t = chrono::high_resolution_clock;

    constexpr size_t frames = 1000;

    const auto seed_val = mt19937::default_seed;
    mt19937 rng( seed_val );

    emitter_t emitter;
    emitter.generate( particle_count, rng );

    auto start = clock_t::now();

    for ( size_t i = 0; i < frames; ++i ) {
        emitter.update();
    }

    auto finish = clock_t::now();
    return finish - start;
}

int main() {
    auto aos_duration = run_test<aos_emitter_t>();
    auto soa_duration = run_test<soa_emitter_t>();
    auto soa_flat_duration = run_test<soa_emitter_flat_t>();
    auto soa_static_duration = run_test<soa_emitter_static_t>();
    auto soa_block_duration = run_test<soa_emitter_block_t>();
  #ifdef HAVE_LIBFLATARRAY
    auto libflatarray_soa_duration = run_test<libflatarray_soa_emitter_t>();
  #endif
  #ifdef HAVE__M128
    auto soa_sse_duration = run_test<soa_emitter_sse_t>();
    auto soa_sse_opt_duration = run_test<soa_emitter_sse_opt_t>();
  #endif

    cout.imbue(std::locale(""));//for thousands separator.
    cout << "particle_count="<< particle_count<< std::endl;
#if 1
    cout << "AoS in " << aos_duration.count() << " seconds" << std::endl;
    cout << "SoA in " << soa_duration.count() << " seconds" << std::endl;
    cout << "SoA flat in " << soa_flat_duration.count() << " seconds" << std::endl;
    cout << "SoA Static in " << soa_static_duration.count() << " seconds" << std::endl;
    cout << "SoA block in " << soa_block_duration.count() << " seconds" << std::endl;
  #ifdef HAVE_LIBFLATARRAY
    cout << "LibFlatArray SoA in " << libflatarray_soa_duration.count() << " seconds" << std::endl;
  #endif
  #ifdef HAVE__M128
    cout << "SoA SSE in " << soa_sse_duration.count() << " seconds" << std::endl;
    cout << "SoA SSE opt in " << soa_sse_opt_duration.count() << " seconds" << std::endl;
  #endif
#endif
    return 0;
}
