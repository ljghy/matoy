#ifndef MATOY_RANGE_HPP_
#define MATOY_RANGE_HPP_

#include <MatrixDefs.h>
#include <cassert>
#include <initializer_list>
#include <vector>
namespace matoy
{
class Range
{
protected:
    index_t m_st;
    index_t m_ed;
    int32_t m_step;
    std::vector<index_t> m_list;
    enum IndexType
    {
        RANGE,
        LIST
    } m_type;

public:
    static const Range all;

    Range() = default;

    Range(index_t st, index_t ed, int32_t step = 1)
        : m_st(st), m_ed(ed), m_step(step), m_type(RANGE)
    {
    }

    Range(const std::initializer_list<index_t> &list)
        : m_st(0), m_ed(0), m_step(0), m_list(list), m_type(LIST)
    {
    }

    template <typename Iterable>
    Range(Iterable list)
        : m_st(0), m_ed(0), m_step(0), m_type(LIST)
    {
        for (index_t i : list)
            m_list.push_back(i);
        assert(!m_list.empty());
    }

    index_t size() const
    {
        if (m_type == RANGE)
            return std::max(0, 1 + ((int32_t)m_ed - (int32_t)m_st) / m_step);
        return m_list.size();
    }

    index_t operator[](index_t i) const
    {
        if (m_type == RANGE)
        {
            index_t ret = m_st + m_step * i;
            return ret;
        }
        assert(i < m_list.size());
        return m_list[i];
    }

    inline bool continuous() const { return (m_type == RANGE && (m_step == 1 || m_st == m_ed)); }
    bool operator==(const Range &rhs) const
    {
        if (m_type == RANGE)
            return (rhs.m_type == RANGE) && (m_st == rhs.m_st) && (m_ed == rhs.m_ed) && (m_step == rhs.m_step);
        return (m_list == rhs.m_list);
    }

    inline index_t st() const { return m_st; }
    inline index_t ed() const { return m_ed; }
};

const Range Range::all(-1, -1);
} // namespace matoy

#endif