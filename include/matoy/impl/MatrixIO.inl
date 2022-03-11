#ifndef MATOY_IMPL_MATRIX_IO_INL_
#define MATOY_IMPL_MATRIX_IO_INL_

#include <iomanip>

namespace matoy
{

template <typename Ty>
inline std::ostream &operator<<(std::ostream &os, const Matrix<Ty> &mat)
{
    os << std::setprecision(OUTPUT_PRECISION) << std::fixed;
    for (index_t r = 0; r < mat.row(); ++r)
    {
        for (index_t c = 0; c < mat.col(); ++c)
        {
            os << std::right << std::setw(OUTPUT_WIDTH) << mat(r, c);
        }
        os << '\n';
    }
    return os;
}

template <typename Ty>
inline std::ostream &operator<<(std::ostream &os, const Matrix<std::complex<Ty>> &mat)
{
    std::ostringstream oss;
    oss << std::setprecision(OUTPUT_PRECISION) << std::fixed;
    for (index_t r = 0; r < mat.row(); ++r)
    {
        for (index_t c = 0; c < mat.col(); ++c)
        {
            oss << std::real(mat(r, c))
                << ' ' << (std::imag(mat(r, c)) > 0 ? '+' : '-') << ' '
                << std::abs(std::imag(mat(r, c))) << 'i';
            os << std::right << std::setw(2 * OUTPUT_WIDTH)
               << oss.str();
            oss.str("");
        }
        os << '\n';
    }
    return os;
}

} // namespace matoy

#endif