#ifndef MATOY_IMPL_MATRIX_OPERATORS_INL_
#define MATOY_IMPL_MATRIX_OPERATORS_INL_

namespace matoy
{
template <typename Ty>
const Matrix<Ty> Matrix<Ty>::_add(const Matrix<Ty> &rhs) const
{
    assert(row() == rhs.row() && col() == rhs.col());
    Matrix<Ty> ret(row(), col());

    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < col(); ++c)
                ret(r, c) = (*this)(r, c) + rhs(r, c);
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < col(); ++c)
            for (index_t r = 0; r < row(); ++r)
                ret(r, c) = (*this)(r, c) + rhs(r, c);
    }
    return ret;
}

template <typename Ty>
const Matrix<Ty> Matrix<Ty>::_neg() const
{
    Matrix<Ty> ret(row(), col());
    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < col(); ++c)
                ret(r, c) = -(*this)(r, c);
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < col(); ++c)
            for (index_t r = 0; r < row(); ++r)
                ret(r, c) = -(*this)(r, c);
    }
    return ret;
}

template <typename Ty>
const Matrix<Ty> Matrix<Ty>::_sub(const Matrix<Ty> &rhs) const
{
    assert(row() == rhs.row() && col() == rhs.col());
    Matrix<Ty> ret(row(), col());
    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < col(); ++c)
                ret(r, c) = (*this)(r, c) - rhs(r, c);
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < col(); ++c)
            for (index_t r = 0; r < row(); ++r)
                ret(r, c) = (*this)(r, c) - rhs(r, c);
    }
    return ret;
}

template <typename Ty>
const Matrix<Ty> Matrix<Ty>::_mul(const Ty &s) const
{
    Matrix<Ty> ret(row(), col());
    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < col(); ++c)
                ret(r, c) = s * (*this)(r, c);
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < col(); ++c)
            for (index_t r = 0; r < row(); ++r)
                ret(r, c) = s * (*this)(r, c);
    }
    return ret;
}

template <typename Ty>
const Matrix<Ty> Matrix<Ty>::_mul(const Matrix<Ty> &rhs) const
{
    assert(col() == rhs.row());
    Matrix<Ty> ret(row(), rhs.col());

    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < rhs.col(); ++c)
            {
                ret(r, c) = (*this)(r, 0) * rhs(0, c);
                for (index_t k = 1; k < col(); ++k)
                    ret(r, c) += (*this)(r, k) * rhs(k, c);
            }
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < rhs.col(); ++c)
            for (index_t r = 0; r < row(); ++r)
            {
                ret(r, c) = (*this)(r, 0) * rhs(0, c);
                for (index_t k = 1; k < col(); ++k)
                    ret(r, c) += (*this)(r, k) * rhs(k, c);
            }
    }
    return ret;
}

template <typename Ty>
const Matrix<Ty> Matrix<Ty>::_div(const Ty &rhs) const
{
    assert(rhs != 0);
    return _mul(1.0 / rhs);
}

template <typename Ty>
Matrix<Ty> &Matrix<Ty>::_addEq(const Matrix<Ty> &rhs)
{
    assert(row() == rhs.row() && col() == rhs.col());
    copy();
    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < col(); ++c)
                (*this)(r, c) += rhs(r, c);
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < col(); ++c)
            for (index_t r = 0; r < row(); ++r)
                (*this)(r, c) += rhs(r, c);
    }
    return *this;
}

template <typename Ty>
Matrix<Ty> &Matrix<Ty>::_subEq(const Matrix<Ty> &rhs)
{
    assert(row() == rhs.row() && col() == rhs.col());
    copy();
    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < col(); ++c)
                (*this)(r, c) -= rhs(r, c);
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < col(); ++c)
            for (index_t r = 0; r < row(); ++r)
                (*this)(r, c) -= rhs(r, c);
    }
    return *this;
}

template <typename Ty>
Matrix<Ty> &Matrix<Ty>::_mulEq(const Ty &rhs)
{
    assert(row() == rhs.row() && col() == rhs.col());
    copy();
    if (row() > col())
    {
#pragma omp parallel for
        for (index_t r = 0; r < row(); ++r)
            for (index_t c = 0; c < col(); ++c)
                (*this)(r, c) *= rhs;
    }
    else
    {
#pragma omp parallel for
        for (index_t c = 0; c < col(); ++c)
            for (index_t r = 0; r < row(); ++r)
                (*this)(r, c) *= rhs;
    }
    return *this;
}
template <typename Ty>
Matrix<Ty> &Matrix<Ty>::_divEq(const Ty &rhs)
{
    assert(rhs != 0);
    return _mulEq(1.0 / rhs);
}

} // namespace matoy
#endif