#ifndef MATOY_IMPL_MATRIX_CONSTRUCTORS_INL_
#define MATOY_IMPL_MATRIX_CONSTRUCTORS_INL_

namespace matoy
{
template <typename Ty>
Matrix<Ty>::Matrix(index_t r0, index_t r, index_t c0, index_t c, std::shared_ptr<MatrixBuffer<Ty>> b, bool t)
    : m_rowOffset(r0), m_row(r), m_colOffset(c0), m_col(c), m_buf(b), m_transpose(t) {}

template <typename Ty>
Matrix<Ty>::Matrix(index_t r, index_t c)
    : m_rowOffset(0), m_row(r), m_colOffset(0), m_col(c),
      m_buf(new MatrixBuffer<Ty>(r, c)), m_transpose(false) {}

template <typename Ty>
Matrix<Ty>::Matrix(const std::initializer_list<const std::initializer_list<Ty>> &elem)
    : m_rowOffset(0), m_row(0), m_colOffset(0), m_col(0),
      m_buf(new MatrixBuffer<Ty>(elem)), m_transpose(false)
{
    m_row = m_buf->row();
    m_col = m_buf->col();
}
} // namespace matoy

#endif