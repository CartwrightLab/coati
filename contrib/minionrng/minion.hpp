/*
MIT License

Copyright (c) 2020 Reed A. Cartwright <reed@cartwrig.ht>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef MINIONRNG_HPP
#define MINIONRNG_HPP

#include <array>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

#if defined(_WIN64) || defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

#if __cpluscplus >= 201103L
#include <chrono>
#include <random>
#else
#include <ctime>
#endif

namespace minion {

class SeedSeq32;
class Random;

namespace detail {

/*
C++11 compatible PRNG engine implementing xoshiro256**

This is xoshiro256** 1.0, our all-purpose, rock-solid generator. It has
excellent (sub-ns) speed, a state (256 bits) that is large enough for
any parallel application, and it passes all tests we are aware of.

The state must be seeded so that it is not everywhere zero. If you have
a 64-bit seed, we suggest to seed a splitmix64 generator and use its
output to fill s.

Based on code written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
*/

class Xoshiro256StarStarEngine {
   public:
    using result_type = uint64_t;
    using state_type = std::array<result_type, 4>;

    explicit Xoshiro256StarStarEngine(const state_type &seed = {0, 0, 0, 0}) { InitState(seed); }

    result_type operator()() { return Next(); }

    void Discard(uint64_t z);

    static constexpr result_type min() { return std::numeric_limits<result_type>::min(); }
    static constexpr result_type max() { return std::numeric_limits<result_type>::max(); }

    friend bool operator==(const Xoshiro256StarStarEngine &left, const Xoshiro256StarStarEngine &right);

    const state_type &state() const { return state_; }
    void SetState(const state_type &state) { state_ = state; };

    void InitState(const state_type &seeds);

   protected:
    result_type Next();

   private:
    state_type state_{};
};

inline uint64_t rotl(const uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

// Advance the engine to the next state and return a pseudo-random value
inline Xoshiro256StarStarEngine::result_type Xoshiro256StarStarEngine::Next() {
    const uint64_t result_starstar = detail::rotl(state_[1] * 5, 7) * 9;

    const uint64_t t = state_[1] << 17;

    state_[2] ^= state_[0];
    state_[3] ^= state_[1];
    state_[1] ^= state_[2];
    state_[0] ^= state_[3];

    state_[2] ^= t;

    state_[3] = detail::rotl(state_[3], 45);

    return result_starstar;
}

// Seed the state of the engine
inline void Xoshiro256StarStarEngine::InitState(const state_type &seeds) {
    // start with well mixed bits
    state_ = {UINT64_C(0x5FAF84EE2AA04CFF), UINT64_C(0xB3A2EF3524D89987), UINT64_C(0x5A82B68EF098F79D),
              UINT64_C(0x5D7AA03298486D6E)};
    // add in the seeds
    state_[0] += seeds[0];
    state_[1] += seeds[1];
    state_[2] += seeds[2];
    state_[3] += seeds[3];
    // check to see if state is all zeros and fix
    if(state_[0] == 0 && state_[1] == 0 && state_[2] == 0 && state_[3] == 0) {
        state_[1] = UINT64_C(0x1615CA18E55EE70C);
    }
    // burn in 256 values
    Discard(256);
}

// Read z values from the Engine and discard
inline void Xoshiro256StarStarEngine::Discard(uint64_t z) {
    for(; z != UINT64_C(0); --z) {
        Next();
    }
}

inline bool operator==(const Xoshiro256StarStarEngine &left, const Xoshiro256StarStarEngine &right) {
    return left.state_ == right.state_;
}

inline bool operator!=(const Xoshiro256StarStarEngine &left, const Xoshiro256StarStarEngine &right) {
    return !(left == right);
}

using RandomEngine = Xoshiro256StarStarEngine;

inline int64_t random_i63(uint64_t u) { return u >> 1; }
inline uint32_t random_u32(uint64_t u) { return u >> 32; }
inline int32_t random_i31(uint64_t u) { return u >> 33; }

// uniformly distributed between [0,max_value)
// Algorithm 5 from Lemire (2018) https://arxiv.org/pdf/1805.10941.pdf
template <typename callback>
uint64_t random_u64_limited(uint64_t max_value, callback &get) {
    uint64_t x = get();
    __uint128_t m = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(max_value);
    auto l = static_cast<uint64_t>(m);
    if(l < max_value) {
        uint64_t t = -max_value % max_value;
        while(l < t) {
            x = get();
            m = static_cast<__uint128_t>(x) * static_cast<__uint128_t>(max_value);
            l = static_cast<uint64_t>(m);
        }
    }
    return m >> 64;
}

inline std::pair<uint32_t, uint32_t> random_u32_pair(uint64_t u) { return {u, u >> 32}; }

// sanity check
static_assert(__FLOAT_WORD_ORDER == __BYTE_ORDER,
              "random_double52 is not implemented if double and uint64_t have different byte orders");

inline double random_f52(uint64_t u) {
    u = (u >> 12) | UINT64_C(0x3FF0000000000000);
    double d;
    std::memcpy(&d, &u, sizeof(d));
    // d - (1.0-(DBL_EPSILON/2.0));
    return d - 0.99999999999999988;
}

inline double random_f53(uint64_t u) {
    auto n = static_cast<int64_t>(u >> 11);
    return n / 9007199254740992.0;
}

inline double random_exp_inv(double f) { return -log(f); }

extern const std::array<int64_t, 256> ek;
extern const std::array<double, 256> ew;
extern const std::array<double, 256> ef;

template <typename callback>
double random_exp_zig_internal(int64_t a, int b, callback &get) {
    constexpr double r = 7.69711747013104972;
    do {
        if(b == 0) {
            return r + random_exp_inv(random_f52(get()));
        }
        double x = a * ew[b];
        if(ef[b - 1] + random_f52(get()) * (ef[b] - ef[b - 1]) < exp(-x)) {
            return x;
        }
        a = random_i63(get());
        b = static_cast<int>(a & 255);
    } while(a > ek[b]);
    return a * ew[b];
}

template <typename callback>
inline double random_exp_zig(callback &get) {
    int64_t a = random_i63(get());
    auto b = static_cast<int>(a & 255);
    if(a <= ek[b]) {
        return a * ew[b];
    }
    return random_exp_zig_internal(a, b, get);
}

}  // namespace detail

// code sanity check
static_assert(std::is_same<uint64_t, detail::RandomEngine::result_type>::value,
              "The result type of RandomEngine is not a uint64_t.");

class Random : public detail::RandomEngine {
   public:
    using detail::RandomEngine::RandomEngine;

    uint64_t bits();
    uint64_t bits(int b);

    uint64_t u64();
    uint64_t u64(uint64_t max_value);

    uint32_t u32();
    std::pair<uint32_t, uint32_t> u32_pair();

    double f52();
    double f53();

    double exp(double mean = 1.0);

    void Seed(const SeedSeq32 &ss);
    void Seed(const uint32_t s);
};

// uniformly distributed between [0,2^64)
inline uint64_t Random::bits() { return detail::RandomEngine::operator()(); }
// uniformly distributed between [0,2^b)
inline uint64_t Random::bits(int b) { return bits() >> (64 - b); }

// uniformly distributed between [0,2^64)
inline uint64_t Random::u64() { return bits(); }

// uniformly distributed between [0,max_value)
inline uint64_t Random::u64(uint64_t max_value) { return detail::random_u64_limited(max_value, *this); }

// uniformly distributed between [0,2^32)
inline uint32_t Random::u32() { return detail::random_u32(bits()); }

// uniformly distributed pair between [0,2^32)
inline std::pair<uint32_t, uint32_t> Random::u32_pair() { return detail::random_u32_pair(bits()); }

// uniformly distributed between (0,1.0)
inline double Random::f52() { return detail::random_f52(bits()); }

// uniformly distributed between [0,1.0)
inline double Random::f53() { return detail::random_f53(bits()); }

// exponential random value with specified mean. mean=1.0/rate
inline double Random::exp(double mean) { return detail::random_exp_zig(*this) * mean; }

namespace detail {

/*
This is a fixed-increment version of Java 8's SplittableRandom generator
See http://dx.doi.org/10.1145/2714064.2660195 and
http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

Based on code written in 2015 by Sebastiano Vigna (vigna@acm.org)
*/

inline uint64_t splitmix64(uint64_t *state) {
    assert(state != nullptr);
    *state += UINT64_C(0x9e3779b97f4a7c15);
    uint64_t z = *state;
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

}  // namespace detail

// Create a sequence of 32-bit states based on hashing 32-bit seeds
class SeedSeq32 {
   public:
    // Stores a vector of seeds in the object
    template <typename... Args>
    explicit SeedSeq32(Args &&... args) : seq_(std::forward<Args>(args)...) {}

    // Generates a range of 32-bit states based on the stored seeds
    // Uses a unique multilinear hash per state (https://arxiv.org/pdf/1202.4961.pdf)
    // Hash is based on a sequence of random 64-bit numbers genrate by splitmix64.
    // splitmix64 is seeded by the initial values in the range
    //
    // Multilinear hash is (m_0 + sum(m_i*u_i) mod 2^64) / 2^32
    // m = buffer of 64-bit unsigned random values
    // u = 32-bit input values that are being hashed

    template <typename It1, typename It2>
    void Generate(It1 first, It2 last) const {
        for(It1 it = first; it != last; ++it) {
            auto s = static_cast<uint64_t>(*it);
            uint64_t sum = detail::splitmix64(&s);
            for(uint32_t u : seq_) {
                sum += detail::splitmix64(&s) * u;
            }
            // If seq_ ends in a zero, the hash is not unique.
            // Add a final value to ensure that this doesn't happen.
            sum += detail::splitmix64(&s) * 1;

            // final value
            *it = static_cast<uint32_t>(sum >> 32);
        }
    }

    // Fills a range of states based on the stored seeds
    template <typename Range>
    void Generate(Range *range) const {
        assert(range != nullptr);
        Generate(std::begin(*range), std::end(*range));
    }

    // Generates a uint32_t seed based on the stored seeds
    uint32_t GenerateOne(uint32_t s = 0xFD57D105u) const {
        std::array<uint32_t, 1> u{s};
        Generate(&u);
        return u[0];
    }

   private:
    std::vector<uint32_t> seq_;
};

inline void Random::Seed(uint32_t s) {
    SeedSeq32 ss(1, s);
    Seed(ss);
}

inline void Random::Seed(const SeedSeq32 &ss) {
    // checking for code assumptions
    static_assert(sizeof(state_type::value_type) == 8, "state_type::value_type does not hold 64 bits");
    static_assert(state_type{}.size() == 4, "state_type does not contain 4 values");

    // starting values / hash seeds
    std::array<uint32_t, 8> values = {0x9272B87Fu, 0xD9F64D09u, 0x6640D56Cu, 0x8CDA60ACu,
                                      0xDEED25EDu, 0x8495FC63u, 0xAEA86A02u, 0x9F129AB9u};
    ss.Generate(values.begin(), values.end());

    // copy 8 32-bit seeds to 4 64-bit seeds
    state_type seeds;
    std::memcpy(seeds.data(), values.data(), seeds.size() * sizeof(seeds[0]));

    InitState(seeds);
}

inline SeedSeq32 create_seed_seq() {
    std::vector<uint32_t> ret;

    // 1. current time
    // 2. current pid
    // 3. 64 random bits (if properly implemented)
    // 4. 64 well-mixed bits on the end (just in case)

#if __cpluscplus >= 201103L
    uint64_t u = std::chrono::high_resolution_clock::now().time_since_epoch().count();
#else
    uint64_t u = time(nullptr);
#endif
    ret.push_back(static_cast<uint32_t>(u));
    ret.push_back(static_cast<uint32_t>(u >> 32));

#if defined(_WIN64) || defined(_WIN32)
    ret.push_back(_getpid());
#else
    ret.push_back(getpid());
#endif

#if __cpluscplus >= 201103L
    ret.push_back(std::random_device{}());
    ret.push_back(std::random_device{}());
#endif

    ret.push_back(0xC8F978DBu);
    ret.push_back(0x0B32F62Eu);

    return SeedSeq32{std::move(ret)};
}

class AliasTable {
   public:
    AliasTable() = default;

    template <typename... Args>
    explicit AliasTable(Args &&... args) {
        create(std::forward<Args>(args)...);
    }

    // create the alias table
    void CreateInplace(std::vector<double> *v);

    // create the alias table
    template <typename... Args>
    void Create(Args &&... args) {
        std::vector<double> vv(std::forward<Args>(args)...);
        CreateInplace(&vv);
    }

    int Get(uint64_t u) const {
        auto yx = detail::random_u32_pair(u >> shr_);
        return (yx.first < p_[yx.second]) ? yx.second : a_[yx.second];
    }

    const std::vector<uint32_t> &a() const { return a_; }
    const std::vector<uint32_t> &p() const { return p_; }

    int operator()(uint64_t u) const { return Get(u); }

   private:
    template <typename T>
    inline static std::pair<T, int> RoundUp(T x) {
        T y = static_cast<T>(2);
        int k = 1;
        for(; y < x; y *= 2, ++k) {
            /*noop*/;
        }
        return std::make_pair(y, k);
    }

    int shr_{0};
    std::vector<uint32_t> a_;
    std::vector<uint32_t> p_;
};

inline void AliasTable::CreateInplace(std::vector<double> *v) {
    assert(v != nullptr);
    assert(v->size() <= std::numeric_limits<uint32_t>::max());
    // round the size of vector up to the nearest power of two
    auto ru = RoundUp(v->size());
    size_t sz = ru.first;
    v->resize(sz, 0.0);
    a_.resize(sz, 0);
    p_.resize(sz, 0);
    // use the number of bits to calculate the right shift operand
    shr_ = 64 - ru.second;

    // find scale for input vector
    double d = std::accumulate(v->begin(), v->end(), 0.0) / sz;

    // find first large and small values
    //     g: current large value index
    //     m: current small value index
    //    mm: next possible small value index
    size_t g, m, mm;
    for(g = 0; g < sz && (*v)[g] < d; ++g) {
        /*noop*/;
    }
    for(m = 0; m < sz && (*v)[m] >= d; ++m) {
        /*noop*/;
    }
    mm = m + 1;

    // construct table
    while(g < sz && m < sz) {
        assert((*v)[m] < d);
        p_[m] = static_cast<uint32_t>(4294967296.0 / d * (*v)[m]);
        a_[m] = static_cast<uint32_t>(g);
        (*v)[g] = ((*v)[g] + (*v)[m]) - d;
        if((*v)[g] >= d || mm <= g) {
            for(m = mm; m < sz && (*v)[m] >= d; ++m) {
                /*noop*/;
            }
            mm = m + 1;
        } else {
            m = g;
        }
        for(; g < sz && (*v)[g] < d; ++g) {
            /*noop*/;
        }
    }
    // if we stopped early fill in the rest
    if(g < sz) {
        p_[g] = std::numeric_limits<uint32_t>::max();
        a_[g] = static_cast<uint32_t>(g);
        for(g = g + 1; g < sz; ++g) {
            if((*v)[g] < d) continue;
            p_[g] = std::numeric_limits<uint32_t>::max();
            a_[g] = static_cast<uint32_t>(g);
        }
    }
    // if we stopped early fill in the rest
    if(m < sz) {
        p_[m] = std::numeric_limits<uint32_t>::max();
        a_[m] = static_cast<uint32_t>(m);
        for(m = mm; m < sz; ++m) {
            if((*v)[m] > d) continue;
            p_[m] = std::numeric_limits<uint32_t>::max();
            a_[m] = static_cast<uint32_t>(m);
        }
    }
}

}  // namespace minion

// MINIONRNG_HPP
#endif
