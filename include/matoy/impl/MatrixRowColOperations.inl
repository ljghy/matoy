#ifndef MATOY_IMPL_MATRIX_ROW_COL_OPERATIONS_INL_
#define MATOY_IMPL_MATRIX_ROW_COL_OPERATIONS_INL_

namespace matoy
{

template <typename Ty>
void Matrix<Ty>::swapRow(index_t r1, index_t r2, index_t st, index_t ed)
{
    if (r1 != r2)
    {
        copy();
        if (ed == -1)
            ed = col();
#pragma omp parallel for
        for (index_t c = st; c < ed; ++c)
            std::swap((*this)(r1, c), (*this)(r2, c));
    }
}

template <typename Ty>
void Matrix<Ty>::swapCol(index_t c1, index_t c2, index_t st, index_t ed)
{
    if (c1 != c2)
    {
        copy();
        if (ed == -1)
            ed = row();
#pragma omp parallel for
        for (index_t r = st; r < ed; ++r)
            std::swap((*this)(r, c1), (*this)(r, c2));
    }
}

template <typename Ty>
void Matrix<Ty>::addRow(index_t src, index_t tar, const Ty &scale, index_t st, index_t ed)
{
    copy();
    if (ed == -1)
        ed = col();
#pragma omp parallel for
    for (index_t c = st; c < ed; ++c)
        (*this)(tar, c) += scale * (*this)(src, c);
}
template <typename Ty>
void Matrix<Ty>::addCol(index_t src, index_t tar, const Ty &scale, index_t st, index_t ed)
{
    copy();
    if (ed == -1)
        ed = row();
#pragma omp parallel for
    for (index_t r = st; r < ed; ++r)
        (*this)(r, tar) += scale * (*this)(r, src);
}

template <typename Ty>
void Matrix<Ty>::mulRow(index_t r, const Ty &k, index_t st, index_t ed)
{
    copy();
    if (ed == -1)
        ed = col();
#pragma omp parallel for
    for (index_t c = st; c < ed; ++c)
        (*this)(r, c) *= k;
}
template <typename Ty>
void Matrix<Ty>::mulCol(index_t c, const Ty &k, index_t st, index_t ed)
{
    copy();
    if (ed == -1)
        ed = row();
#pragma omp parallel for
    for (index_t r = 0; r < ed; ++r)
        (*this)(r, c) *= k;
}

template <typename Ty>
inline Matrix<Ty> rowCat(const Matrix<Ty> &left, const Matrix<Ty> &right)
{
    index_t row = left.row();
    assert(row = right.row());
    index_t coll = left.col(), colr = right.col();
    Matrix<Ty> ret(row, coll + colr);
#pragma omp parallel for
    for (index_t r = 0; r < row; ++r)
    {
        for (index_t c = 0; c < coll; ++c)
            ret(r, c) = left(r, c);
        for (index_t c = 0; c < colr; ++c)
            ret(r, c + coll) = right(r, c);
    }
    return ret;
}

template <typename Ty>
inline Matrix<Ty> colCat(const Matrix<Ty> &up, const Matrix<Ty> &down)
{
    index_t col = up.col();
    assert(col = down.col());
    index_t rowu = up.row(), rowd = down.row();
    Matrix<Ty> ret(rowu + rowd, col);
#pragma omp parallel for
    for (index_t c = 0; c < col; ++c)
    {
        for (index_t r = 0; r < rowu; ++r)
            ret(r, c) = up(r, c);
        for (index_t r = 0; r < rowd; ++r)
            ret(r + rowu, c) = down(r, c);
    }
    return ret;
}

} // namespace matoy

#endif