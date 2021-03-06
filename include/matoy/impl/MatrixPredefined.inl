#ifndef MATOY_IMPL_MATRIX_PREDEFINED_INL_
#define MATOY_IMPL_MATRIX_PREDEFINED_INL_

#include <functional>

namespace matoy
{
template <typename Ty>
inline Matrix<Ty> eye(index_t row, index_t col = -1)
{
    if (col == -1)
        col = row;
    Matrix<Ty> ret(row, col);
    index_t n = std::min(row, col);
    for (index_t r = 0; r < n; ++r)
        for (index_t c = 0; c < n; ++c)
            ret(r, c) = Ty(r == c);
    return ret;
}

template <typename Ty>
inline Matrix<Ty> scalar(const Ty &diag, index_t row, index_t col = -1)
{
    if (col == -1)
        col = row;
    Matrix<Ty> ret(row, col);
    index_t n = std::min(row, col);
    for (index_t r = 0; r < n; ++r)
        for (index_t c = 0; c < n; ++c)
            ret(r, c) = (r == c) ? diag : Ty(0.0);
    return ret;
}

template <typename Ty>
inline Matrix<Ty> rand(const std::function<Ty()> &randGen, index_t row, index_t col = -1)
{
    if (col == -1)
        col = row;
    Matrix<Ty> ret(row, col);
    for (index_t r = 0; r < row; ++r)
        for (index_t c = 0; c < col; ++c)
            ret(r, c) = randGen();
    return ret;
}

template <typename Ty>
inline Matrix<Ty> zeros(index_t row, index_t col = -1)
{
    if (col == -1)
        col = row;
    Matrix<Ty> ret(row, col);
    for (index_t r = 0; r < row; ++r)
        for (index_t c = 0; c < col; ++c)
            ret(r, c) = 0;
    return ret;
}
template <typename Ty>
inline Matrix<Ty> ones(index_t row, index_t col = -1)
{
    if (col == -1)
        col = row;
    Matrix<Ty> ret(row, col);
    for (index_t r = 0; r < row; ++r)
        for (index_t c = 0; c < col; ++c)
            ret(r, c) = 1.0;
    return ret;
}

template <typename Ty>
inline Matrix<Ty> diag(const Matrix<Ty> &vec, index_t row, index_t col = -1)
{
    assert(vec.col() == 1);
    Matrix<Ty> ret(zeros<Ty>(row, col));
    index_t n = std::min(std::min(vec.row(), row), col);
    for (index_t i = 0; i < n; ++i)
        ret(i, i) = vec(i, 0);
    return ret;
}

} // namespace matoy

#endif