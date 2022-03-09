#ifndef MATOY_IMPL_MATRIX_DET_INL_
#define MATOY_IMPL_MATRIX_DET_INL_

namespace matoy
{
template <typename Ty>
inline Ty det(Matrix<Ty> mat)
{
    assert(mat.isSqr());

    index_t n = mat.row();
    Ty max, d = 1;
    index_t maxRow;
    for (index_t i = 0; i < n; ++i)
    {
        max = std::abs(mat(i, i));
        maxRow = i;
        for (index_t j = i + 1; j < n; ++j)
            if (std::abs(mat(j, i)) > max)
            {
                max = std::abs(mat(j, i));
                maxRow = j;
            }
        if (IS_ZERP(max))
            return 0;
        if (i != maxRow)
        {
            mat.swapRow(i, maxRow, i);
            d = -d;
        }
        for (index_t j = i + 1; j < n; ++j)
            mat.addRow(i, j, -mat(j, i) / mat(i, i), i);
    }
    for (index_t i = 0; i < n; ++i)
        d *= mat(i, i);
    return d;
}
} // namespace matoy

#endif