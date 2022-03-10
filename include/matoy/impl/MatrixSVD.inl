#ifndef MATOY_IMPL_MATRIX_SVD_INL_
#define MATOY_IMPL_MATRIX_SVD_INL_

namespace matoy
{

template <typename Ty>
inline std::tuple<Matrix<Ty>, Matrix<Ty>, Matrix<Ty>>
svd(const Matrix<Ty> &mat, double tol = DEFAULT_TOL, index_t iterCount = DEFAULT_ITER_COUNT)
{
    index_t row = mat.row(), col = mat.col(), n = row + col;
    Matrix<Ty> b = zeros<Ty>(n);
    b.setSubmat({col, n - 1}, {0, col - 1}, mat);
    b.setSubmat({0, col - 1}, {col, n - 1}, mat.T());

    auto p = eig(b, tol, iterCount);

    // TODO
}
} // namespace matoy

#endif