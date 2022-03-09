#ifndef MATOY_MATRIX_BUFFER_HPP_
#define MATOY_MATRIX_BUFFER_HPP_

#include <MatrixDefs.h>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <vector>
namespace matoy
{
template <typename Ty>
class MatrixBuffer
{
protected:
    index_t m_row;
    index_t m_col;
    std::vector<Ty> m_elem;

public:
    MatrixBuffer(index_t r, index_t c) : m_row(r), m_col(c) { m_elem.resize(r * c); }
    MatrixBuffer(const std::initializer_list<const std::initializer_list<Ty>> &elem)
        : m_row(0), m_col(0)
    {
        assert(elem.size() && elem.begin()->size());
        m_row = elem.size();
        m_col = elem.begin()->size();
        m_elem.reserve(m_row * m_col);
        for (const auto &l : elem)
        {
            assert(l.size() == m_col);
            for (const auto &e : l)
                m_elem.push_back(e);
        }
    }

    MatrixBuffer(const MatrixBuffer<Ty> &) = default;
    MatrixBuffer &operator=(const MatrixBuffer<Ty> &) = default;

    inline index_t row() const { return m_row; }
    inline index_t col() const { return m_col; }

    inline Ty &getElemRef(index_t r0, index_t r, index_t c0, index_t c, bool transpose)
    {
        index_t r1, c1;
        if (transpose)
        {
            r1 = r0 + c;
            c1 = c0 + r;
        }
        else
        {
            r1 = r0 + r;
            c1 = c0 + c;
        }
        assert(r1 < m_row && c1 < m_col);
        return m_elem[m_col * r1 + c1];
    }
    inline const Ty getElemVal(index_t r0, index_t r, index_t c0, index_t c, bool transpose) const
    {
        index_t r1, c1;
        if (transpose)
        {
            r1 = r0 + c;
            c1 = c0 + r;
        }
        else
        {
            r1 = r0 + r;
            c1 = c0 + c;
        }
        assert(r1 < m_row && c1 < m_col);
        return m_elem[m_col * r1 + c1];
    }

    virtual ~MatrixBuffer() {}
};
} // namespace matoy
#endif
