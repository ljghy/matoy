#ifndef MATOY_IMPL_MATRIX_TRANSPOSE_INL_
#define MATOY_IMPL_MATRIX_TRANSPOSE_INL_

namespace matoy
{
template <typename Ty>
Matrix<Ty> Matrix<Ty>::T() const
{
    auto ret = *this;
    ret.transpose();
    return ret;
}

template <typename Ty>
void Matrix<Ty>::transpose()
{
    m_transpose = !m_transpose;
}

template <>
void Matrix<std::complex<float>>::transpose()
{
    m_transpose = !m_transpose;
    for (index_t r = 0; r < row(); ++r)
        for (index_t c = 0; c < col(); ++c)
            (*this)(r, c) = std::conj((*this)(r, c));
}

template <>
void Matrix<std::complex<double>>::transpose()
{
    m_transpose = !m_transpose;
    for (index_t r = 0; r < row(); ++r)
        for (index_t c = 0; c < col(); ++c)
            (*this)(r, c) = std::conj((*this)(r, c));
}

template <>
void Matrix<std::complex<long double>>::transpose()
{
    m_transpose = !m_transpose;
    for (index_t r = 0; r < row(); ++r)
        for (index_t c = 0; c < col(); ++c)
            (*this)(r, c) = std::conj((*this)(r, c));
}

} // namespace matoy

#endif