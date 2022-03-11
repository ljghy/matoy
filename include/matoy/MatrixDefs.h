#ifndef MATOY_MATRIX_DEFS_H_
#define MATOY_MATRIX_DEFS_H_
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdlib>
namespace matoy
{
#ifndef _MSC_VER
using index_t = uint32_t;
#else
using index_t = int32_t;
#endif
constexpr double DEFAULT_TOL = 1e-14;
constexpr index_t DEFAULT_ITER_COUNT = 64;

constexpr unsigned OUTPUT_PRECISION = 4;
constexpr unsigned OUTPUT_WIDTH = 10;

template <typename Ty>
inline bool IS_ZERO(const Ty &x) { return std::abs(x) < DEFAULT_TOL; }

template <typename Ty>
inline Ty randnum() { return (Ty)std::rand() / RAND_MAX; }

template <>
inline std::complex<float> randnum<std::complex<float>>()
{
    return {(float)std::rand() / RAND_MAX, (float)std::rand() / RAND_MAX};
}

template <>
inline std::complex<double> randnum<std::complex<double>>()
{
    return {(double)std::rand() / RAND_MAX, (double)std::rand() / RAND_MAX};
}

template <>
inline std::complex<long double> randnum<std::complex<long double>>()
{
    return {(long double)std::rand() / RAND_MAX, (long double)std::rand() / RAND_MAX};
}

} // namespace matoy

#endif