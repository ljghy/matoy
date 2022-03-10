#ifndef MATOY_RANGE_HPP_
#define MATOY_RANGE_HPP_

#include <MatrixDefs.h>
#include <cassert>
#include <initializer_list>
namespace matoy
{
class Range
{
protected:
    index_t m_st;
    index_t m_ed;
    int32_t m_step;

public:
    static const Range all;

    Range() = default;

    Range(index_t st, index_t ed, int32_t step = 1)
        : m_st(st), m_ed(ed), m_step(step)
    {
    }

    index_t size() const
    {
        return std::max(0, 1 + ((int32_t)m_ed - (int32_t)m_st) / m_step);
    }

    index_t operator[](index_t i) const
    {
        return m_st + m_step * i;
    }

    inline bool continuous() const { return (m_step == 1 || m_st == m_ed); }
    bool operator==(const Range &rhs) const
    {
        return (m_st == rhs.m_st) && (m_ed == rhs.m_ed) && (m_step == rhs.m_step);
    }

    inline index_t st() const { return m_st; }
    inline index_t ed() const { return m_ed; }
};

const Range Range::all(-1, -1);
} // namespace matoy

#endif