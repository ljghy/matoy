#ifndef MATOY_IMPL_MATRIX_SVD_INL_
#define MATOY_IMPL_MATRIX_SVD_INL_

namespace matoy
{

template <typename Ty>
inline std::tuple<Matrix<Ty>, Matrix<Ty>, Matrix<Ty>>
svd(const Matrix<Ty> &mat, double tol = DEFAULT_TOL, index_t iterCount = DEFAULT_ITER_COUNT)
{
    index_t row = mat.row(), col = mat.col(), n = row + col;

    Matrix<std::complex<Ty>> lam, v;
    std::tie(lam, v) = eig(mat.T() * mat, tol, iterCount);
    Matrix<Ty> lamr = real(lam), vr = real(v);

    for (index_t i = 0; i < col; ++i)
        lamr(i, 0) = std::sqrt(lamr(i, 0));

    for (index_t i = 0; i < col - 1; ++i)
    {
        bool swapped = false;
        for (index_t j = 0; j < col - 1 - i; ++j)
            if (lamr(j, 0) < lamr(j + 1, 0))
            {
                std::swap(lamr(j, 0), lamr(j + 1, 0));
                vr.swapCol(j, j + 1);
                swapped = true;
            }
        if (!swapped)
            break;
    }
    Matrix<Ty> s = diag(lamr, row, col);
    Matrix<Ty> u = (mat * vr)(Range::all, {0, row - 1});
    for (index_t i = 0; i < row; ++i)
        if (!IS_ZERO(lamr(i, 0)))
            u.mulCol(i, (Ty)1.0 / lamr(i, 0));
    return {u, s, vr};
}

template <typename Ty>
inline std::tuple<Matrix<std::complex<Ty>>, Matrix<std::complex<Ty>>, Matrix<std::complex<Ty>>>
svd(const Matrix<std::complex<Ty>> &mat, double tol = DEFAULT_TOL, index_t iterCount = DEFAULT_ITER_COUNT)
{
    index_t row = mat.row(), col = mat.col(), n = row + col;

    Matrix<std::complex<Ty>> lam, v;
    std::tie(lam, v) = eig(mat.T() * mat, tol, iterCount);

    auto lamr = real(lam);
    for (index_t i = 0; i < col; ++i)
        lamr(i, 0) = std::sqrt(lamr(i, 0));

    for (index_t i = 0; i < col - 1; ++i)
    {
        bool swapped = false;
        for (index_t j = 0; j < col - 1 - i; ++j)
            if (lamr(j, 0) < lamr(j + 1, 0))
            {
                std::swap(lamr(j, 0), lamr(j + 1, 0));
                v.swapCol(j, j + 1);
                swapped = true;
            }
        if (!swapped)
            break;
    }
    Matrix<std::complex<Ty>> s = diag(cmplx(lamr), row, col);
    Matrix<std::complex<Ty>> u = (mat * v)(Range::all, {0, row - 1});
    for (index_t i = 0; i < row; ++i)
        if (!IS_ZERO(lamr(i, 0)))
            u.mulCol(i, std::complex<Ty>((Ty)1.0 / lamr(i, 0)));
    return {u, s, v};
}

} // namespace matoy

#endif