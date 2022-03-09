#ifndef MATOY_IMPL_MATRIX_QR_INL_
#define MATOY_IMPL_MATRIX_QR_INL_
#include <cstdio>
#include <utility>
namespace matoy
{

template <typename Ty>
inline Matrix<Ty> householder(const Matrix<Ty> &v)
{
    assert(v.col() == 1);
    return eye<Ty>(v.row()) - 2.0 / sqrNorm(v) * v * v.T();
}

template <typename Ty>
inline std::array<Matrix<Ty>, 2> qr(const Matrix<Ty> &mat)
{
    index_t row = mat.row(), col = mat.col();
    Matrix<Ty> I(eye<Ty>(row));
    Matrix<Ty> Q(I), H(I), R(mat);

    for (index_t i = 0; i < std::min(row, col); ++i)
    {
        Matrix<Ty> x = R(Range(i, row - 1), i);
        Matrix<Ty> w = orthoBasis<Ty>(row - i, 0, norm(x));
        if (!IS_ZERO(x(0, 0)))
            w(0, 0) *= x(0, 0) / std::abs(x(0, 0));
        Matrix<Ty> v1 = w - x, v2 = w + x;
        auto s1 = sqrNorm(v1), s2 = sqrNorm(v2);
        if (IS_ZERO(std::max(s1, s2)))
            continue;
        H.setSubmat(Range(i, row - 1), Range(i, row - 1), householder(s1 > s2 ? v1 : v2));
        Q = Q * H;
        R = H * R;
        R.setSubmat(Range(i + 1, row - 1), i, zeros<Ty>(row - i - 1, 1));
        H = I;
    }
    return {Q, R};
}

} // namespace matoy
#endif