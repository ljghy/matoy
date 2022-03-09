#ifndef MATOY_IMPL_MATRIX_INDEX_INL_
#define MATOY_IMPL_MATRIX_INDEX_INL_
namespace matoy
{
template <typename Ty>
Ty &Matrix<Ty>::operator()(index_t r, index_t c)
{
    copy();
    return m_buf->getElemRef(m_rowOffset, r, m_colOffset, c, m_transpose);
}
template <typename Ty>
const Ty Matrix<Ty>::operator()(index_t r, index_t c) const
{
    return m_buf->getElemVal(m_rowOffset, r, m_colOffset, c, m_transpose);
}
template <typename Ty>
const Matrix<Ty> Matrix<Ty>::operator()(const Range &rRange, const Range &cRange) const
{
    return _getByIndexRange(rRange, cRange);
}
template <typename Ty>
const Matrix<Ty> Matrix<Ty>::operator()(index_t r, const Range &cRange) const
{
    return _getByIndexRange(Range(r, r), cRange);
}
template <typename Ty>
const Matrix<Ty> Matrix<Ty>::operator()(const Range &rRange, index_t c) const
{
    return _getByIndexRange(rRange, Range(c, c));
}
template <typename Ty>
void Matrix<Ty>::setSubmat(const Range &rRange, const Range &cRange, const Matrix<Ty> &m)
{
    _setByIndexRange(rRange, cRange, m);
}
template <typename Ty>
void Matrix<Ty>::setSubmat(const index_t &r, const Range &cRange, const Matrix<Ty> &m)
{
    _setByIndexRange(Range(r, r), cRange, m);
}
template <typename Ty>
void Matrix<Ty>::setSubmat(const Range &rRange, const index_t &c, const Matrix<Ty> &m)
{
    _setByIndexRange(rRange, Range(c, c), m);
}
template <typename Ty>
const Matrix<Ty> Matrix<Ty>::_getByIndexRange(Range rowRange, Range colRange) const
{
    if (rowRange == Range::all)
        rowRange = Range(0, row() - 1);
    if (colRange == Range::all)
        colRange = Range(0, col() - 1);
    if (rowRange.continuous() && colRange.continuous())

        return m_transpose ? Matrix<Ty>{m_rowOffset + colRange.st(), colRange.size(),
                                        m_colOffset + rowRange.st(), rowRange.size(),
                                        m_buf, m_transpose}
                           : Matrix<Ty>{m_rowOffset + rowRange.st(), rowRange.size(),
                                        m_colOffset + colRange.st(), colRange.size(),
                                        m_buf, m_transpose};

    Matrix<Ty> ret(rowRange.size(), colRange.size());
#pragma omp parallel for
    for (index_t r = 0; r < rowRange.size(); ++r)
        for (index_t c = 0; c < colRange.size(); ++c)
            ret(r, c) = (*this)(rowRange[r], colRange[c]);
    return ret;
}
template <typename Ty>
void Matrix<Ty>::_setByIndexRange(Range rowRange, Range colRange, const Matrix<Ty> &m)
{
    assert(rowRange.size() == m.row() && colRange.size() == m.col());
    copy();
    if (rowRange == Range::all)
        rowRange = Range(0, row() - 1);
    if (colRange == Range::all)
        colRange = Range(0, col() - 1);
#pragma omp parallel for
    for (index_t r = 0; r < rowRange.size(); ++r)
        for (index_t c = 0; c < colRange.size(); ++c)
            (*this)(rowRange[r], colRange[c]) = m(r, c);
}
} // namespace matoy

#endif