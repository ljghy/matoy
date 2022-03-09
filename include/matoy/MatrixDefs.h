#ifndef MATOY_MATRIX_DEFS_H_
#define MATOY_MATRIX_DEFS_H_
#include <cmath>
#include <complex>
#include <cstdint>
namespace matoy
{
#ifndef _MSC_VER
using index_t = uint32_t;
#else
using index_t = int32_t;
#endif
constexpr double DEFAULT_TOL = 1e-14;

template <typename Ty>
inline bool IS_ZERO(const Ty &x) { return std::abs(x) < DEFAULT_TOL; }
} // namespace matoy

#endif