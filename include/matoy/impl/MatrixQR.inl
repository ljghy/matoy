#ifndef MATOY_IMPL_MATRIX_QR_INL_
#define MATOY_IMPL_MATRIX_QR_INL_
#include <cstdio>
namespace matoy
{

template <typename Ty>
inline Matrix<Ty> householder(const Matrix<Ty> &v)
{
    assert(v.col() == 1);
    return eye<Ty>(v.row()) - 2.0 / sqrNorm(v) * v * v.T();
}

template <typename Ty>
inline std::tuple<Matrix<Ty>, Matrix<Ty>> qr(const Matrix<Ty> &mat)
{
    index_t row = mat.row(), col = mat.col();
    Matrix<Ty> q(eye<Ty>(row)), r(mat);

    for (index_t i = 0; i < std::min(row, col); ++i)
    {
        Matrix<Ty> x = r({i, row - 1}, i);
        Matrix<Ty> w = orthoBasis<Ty>(row - i, 0, norm(x));
        if (!IS_ZERO(x(0, 0)))
            w(0, 0) *= x(0, 0) / std::abs(x(0, 0));
        Matrix<Ty> v1 = w - x, v2 = w + x;
        auto s1 = sqrNorm(v1), s2 = sqrNorm(v2);
        if (IS_ZERO(std::max(s1, s2)))
            continue;

        if (s1 < s2)
        {
            w = -w;
            v1 = v2;
        }
        auto h = householder(v1);
        if (i > 0)
            q.setSubmat({0, i - 1}, {i, row - 1}, q({0, i - 1}, {i, row - 1}) * h);
        q.setSubmat({i, row - 1}, {i, row - 1}, q({i, row - 1}, {i, row - 1}) * h);

        r.setSubmat({i, row - 1}, i, w);
        r.setSubmat({i, row - 1}, {i + 1, col - 1}, h * r({i, row - 1}, {i + 1, col - 1}));
    }
    return {q, r};
}

} // namespace matoy
#endif