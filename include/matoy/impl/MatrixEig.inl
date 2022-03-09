#ifndef MATOY_IMPL_MATRIX_EIG_INL_
#define MATOY_IMPL_MATRIX_EIG_INL_

namespace matoy
{

template <typename Ty>
inline std::pair<Ty, Matrix<Ty>> powit(const Matrix<Ty> &mat, index_t iterCount = DEFAULT_ITER_COUNT)
{
    assert(mat.isSqr());
    Matrix<Ty> x = rand<Ty>(randnum<Ty>, mat.row(), 1), u(x);
    Ty lam;
    for (index_t i = 0; i < iterCount; ++i)
    {
        u = normalized(x);
        x = mat * u;
        lam = (u.T() * x).toNum();
    }
    return {lam, normalized(u)};
}

} // namespace matoy
#endif