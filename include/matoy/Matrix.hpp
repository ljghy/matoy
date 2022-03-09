#ifndef MATOY_MATRIX_HPP_
#define MATOY_MATRIX_HPP_

#include <MatrixBuffer.hpp>
#include <Range.hpp>
#include <array>
#include <cmath>
#include <initializer_list>
#include <memory>
#include <omp.h>

namespace matoy
{

template <typename Ty>
class Matrix
{
protected:
    index_t m_rowOffset;
    index_t m_row;
    index_t m_colOffset;
    index_t m_col;
    std::shared_ptr<MatrixBuffer<Ty>> m_buf;
    bool m_transpose;

public:
    Matrix<Ty> T() const;
    void transpose();

    void swapRow(index_t r1, index_t r2, index_t st = 0, index_t ed = -1);
    void swapCol(index_t c1, index_t c2, index_t st = 0, index_t ed = -1);
    void addRow(index_t src, index_t tar, const Ty &scale, index_t st = 0, index_t ed = -1);
    void addCol(index_t src, index_t tar, const Ty &scale, index_t st = 0, index_t ed = -1);
    void mulRow(index_t r, const Ty &k, index_t st = 0, index_t ed = -1);
    void mulCol(index_t c, const Ty &k, index_t st = 0, index_t ed = -1);

public:
    const Matrix<Ty> operator+(const Matrix<Ty> &rhs) const { return _add(rhs); }
    const Matrix<Ty> operator-() const { return _neg(); }
    const Matrix<Ty> operator-(const Matrix<Ty> &rhs) const { return _sub(rhs); }
    const Matrix<Ty> operator*(const Matrix<Ty> &rhs) const { return _mul(rhs); }
    const Matrix<Ty> operator*(const Ty &rhs) const { return _mul(rhs); }
    friend Matrix<Ty> operator*(const Ty &lhs, const Matrix<Ty> &m) { return m._mul(lhs); }
    const Matrix<Ty> operator/(const Ty &rhs) const { return _div(rhs); }

    inline Ty toNum() const
    {
        assert(row() == 1 && col() == 1);
        return (*this)(0, 0);
    }

    Matrix<Ty> &operator+=(const Matrix<Ty> &rhs) { return _addEq(rhs); }
    Matrix<Ty> &operator-=(const Matrix<Ty> &rhs) { return _subEq(rhs); }
    Matrix<Ty> &operator*=(const Ty &rhs) { return _mulEq(rhs); }
    Matrix<Ty> &operator/=(const Ty &rhs) { return _divEq(rhs); }

protected:
    Matrix(index_t r0, index_t r, index_t c0, index_t c, std::shared_ptr<MatrixBuffer<Ty>> b, bool t);

    inline void copy()
    {
        if (m_buf.use_count() > 1)
            m_buf.reset(new MatrixBuffer<Ty>(*m_buf));
    }

public:
    Matrix(index_t r, index_t c);
    Matrix(const std::initializer_list<const std::initializer_list<Ty>> &elem);

    Matrix(const Matrix<Ty> &) = default;
    Matrix &operator=(const Matrix<Ty> &) = default;

    inline index_t row() const { return m_transpose ? m_col : m_row; }
    inline index_t col() const { return m_transpose ? m_row : m_col; }
    inline bool isSqr() const { return m_row == m_col; }

    Ty &operator()(index_t r, index_t c);
    const Ty operator()(index_t r, index_t c) const;
    const Matrix<Ty> operator()(const Range &rRange, const Range &cRange) const;
    const Matrix<Ty> operator()(index_t r, const Range &cRange) const;
    const Matrix<Ty> operator()(const Range &rRange, index_t c) const;

    void setSubmat(const Range &rRange, const Range &cRange, const Matrix<Ty> &m);
    void setSubmat(const index_t &r, const Range &cRange, const Matrix<Ty> &m);
    void setSubmat(const Range &rRange, const index_t &c, const Matrix<Ty> &m);

    virtual ~Matrix() = default;

protected:
    const Matrix<Ty> _add(const Matrix<Ty> &rhs) const;
    const Matrix<Ty> _neg() const;
    const Matrix<Ty> _sub(const Matrix<Ty> &rhs) const;
    const Matrix<Ty> _mul(const Ty &s) const;
    const Matrix<Ty> _mul(const Matrix<Ty> &rhs) const;
    const Matrix<Ty> _div(const Ty &rhs) const;

    Matrix<Ty> &_addEq(const Matrix<Ty> &rhs);
    Matrix<Ty> &_subEq(const Matrix<Ty> &rhs);
    Matrix<Ty> &_mulEq(const Ty &rhs);
    Matrix<Ty> &_divEq(const Ty &rhs);

    const Matrix<Ty> _getByIndexRange(Range rowRange, Range colRange) const;
    void _setByIndexRange(Range rowRange, Range colRange, const Matrix<Ty> &m);
};

} // namespace matoy

#include <impl/MatrixConstructors.inl>
#include <impl/MatrixDet.inl>
#include <impl/MatrixEig.inl>
#include <impl/MatrixIndex.inl>
#include <impl/Vector.inl>
#include <impl/MatrixOperators.inl>
#include <impl/MatrixPredefined.inl>
#include <impl/MatrixQR.inl>
#include <impl/MatrixRREF.inl>
#include <impl/MatrixRowColOperations.inl>
#include <impl/MatrixTranspose.inl>

#endif