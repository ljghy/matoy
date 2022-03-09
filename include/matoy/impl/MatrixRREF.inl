#ifndef MATOY_IMPL_MATRIX_RREF_INL_
#define MATOY_IMPL_MATRIX_RREF_INL_
namespace matoy
{
template <typename Ty>
inline Matrix<Ty> rref(Matrix<Ty> mat)
{
    index_t row = mat.row(), col = mat.col();
    double max;
    index_t maxRow, r = 0;
    for (index_t c = 0; c < col && r < row; ++c)
    {
        max = std::abs(mat(r, c));
        if (IS_ZERO(max))
            mat(r, c) = max = 0;
        maxRow = r;
        for (index_t i = r + 1; i < row; ++i)
        {
            if (IS_ZERO(mat(i, c)))
                mat(i, c) = 0;
            else if (std::abs(mat(i, c)) > max)
            {
                max = std::abs(mat(i, c));
                maxRow = i;
            }
        }
        if (IS_ZERO(max))
            continue;
        mat.swapRow(r, maxRow, c);
        for (index_t i = r + 1; i < row; ++i)
        {
            mat.addRow(r, i, -mat(i, c) / mat(r, c), c + 1);
            mat(i, c) = 0;
        }
        ++r;
    }

    index_t i = row, c = 0;
    do
    {
        --i;
        while (c < col && IS_ZERO(mat(i, c)))
            ++c;
        if (c == col)
        {
            c = 0;
            continue;
        }
        mat.mulRow(i, 1.0 / mat(i, c), c + 1);
        mat(i, c) = 1;
        for (index_t j = 0; j < i; ++j)
        {
            mat.addRow(i, j, -mat(j, c), c + 1);
            mat(j, c) = 0;
        }
        c = 0;
    } while (i > 0);

    return mat;
}

} // namespace matoy

#endif