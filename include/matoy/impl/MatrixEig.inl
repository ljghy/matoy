#ifndef MATOY_IMPL_MATRIX_EIG_INL_
#define MATOY_IMPL_MATRIX_EIG_INL_

namespace matoy
{

template <typename Ty>
inline Ty powit(const Matrix<Ty> &mat, index_t iterCount = DEFAULT_ITER_COUNT)
{
    assert(mat.isSqr());
    Matrix<Ty> x = rand<Ty>(randnum<Ty>, mat.row(), 1), u(x);
    Ty lam;
    for (index_t i = 0; i < iterCount; ++i)
    {
        u = normalized(x);
        x = mat * u;
        lam = u.T() * x;
    }
    return lam;
}

template <typename Ty>
inline std::tuple<Matrix<Ty>, Matrix<Ty>> upperHessenberg(Matrix<Ty> mat)
{
    assert(mat.isSqr());
    index_t n = mat.row();
    Matrix<Ty> q(eye<Ty>(n));
    for (index_t i = 0; i < n - 2; ++i)
    {
        auto x = mat({i + 1, n - 1}, i);
        auto w = orthoBasis<Ty>(n - i - 1, 0, norm(x));
        if (!IS_ZERO(x(0, 0)))
            w(0, 0) *= x(0, 0) / std::abs(x(0, 0));
        Matrix<Ty> v1 = w - x, v2 = w + x;
        auto s1 = sqrNorm(v1), s2 = sqrNorm(v2);
        if (IS_ZERO(std::max(s1, s2)))
            continue;

        if (s1 < s2)
        {
            w = -w;
            v1 = v2;
        }
        auto h = householder(v1);
        mat.setSubmat({i + 1, n - 1}, i, w);
        mat.setSubmat({i + 1, n - 1}, {i + 1, n - 1}, h * mat({i + 1, n - 1}, {i + 1, n - 1}) * h);
        mat.setSubmat({0, i}, {i + 1, n - 1}, mat({0, i}, {i + 1, n - 1}) * h);

        q.setSubmat(Range::all, {i + 1, n - 1}, q(Range::all, {i + 1, n - 1}) * h);
    }
    return {mat, q};
}

template <typename Ty>
inline std::tuple<Matrix<Ty>, Matrix<Ty>> shiftedQR(Matrix<Ty> mat, double tol, index_t iterCount)
{
    assert(mat.isSqr());
    index_t n = mat.row();
    Matrix<Ty> lam(n, 1), q0(eye<Ty>(n));
    while (n > 1)
    {
        index_t count = 0;
        while (++count < iterCount)
        {
            auto s = scalar(mat(n - 1, n - 1), n);
            Matrix<Ty> q, r;
            std::tie(q, r) = qr(mat - s);
            q0.setSubmat(Range::all, {0, n - 1}, q0(Range::all, {0, n - 1}) * q);
            mat = r * q + s;

            auto max = std::abs(mat(n - 1, 0));
            for (index_t i = 1; i < n - 1; ++i)
                if (std::abs(mat(n - 1, i)) > max)
                    max = std::abs(mat(n - 1, i));
            if (max < tol)
                break;
        }
        if (count < iterCount)
        {
            lam(n - 1, 0) = mat(n - 1, n - 1);
            --n;
            mat = mat({0, n - 1}, {0, n - 1});
        }
        else
        {
            auto tr = mat(n - 2, n - 2) + mat(n - 1, n - 1);
            auto d = std::pow((mat(n - 2, n - 2) - mat(n - 1, n - 1)), 2) +
                     Ty(4.0) * mat(n - 1, n - 2) * mat(n - 2, n - 1);
            lam(n - 2, 0) = (tr + std::sqrt(d)) * Ty(0.5);
            lam(n - 1, 0) = (tr - std::sqrt(d)) * Ty(0.5);
            n -= 2;
            mat = mat({0, n - 1}, {0, n - 1});
        }
    }
    if (n > 0)
        lam(0, 0) = mat(0, 0);
    return {lam, q0};
}

template <typename Ty>
inline std::tuple<Matrix<std::complex<Ty>>, Matrix<std::complex<Ty>>>
eig(Matrix<Ty> mat, double tol = DEFAULT_TOL, index_t iterCount = DEFAULT_ITER_COUNT)
{
    Matrix<Ty> u, q1;
    std::tie(u, q1) = upperHessenberg(mat);
    Matrix<std::complex<Ty>> lam, q2;
    std::tie(lam, q2) = shiftedQR(cmplx(u), tol, iterCount);
    return {lam, cmplx(q1) * q2};
}

template <typename Ty>
inline std::tuple<Matrix<std::complex<Ty>>, Matrix<std::complex<Ty>>>
eig(Matrix<std::complex<Ty>> mat, double tol = DEFAULT_TOL, index_t iterCount = DEFAULT_ITER_COUNT)
{
    Matrix<std::complex<Ty>> u, q1;
    std::tie(u, q1) = upperHessenberg(mat);
    Matrix<std::complex<Ty>> lam, q2;
    std::tie(lam, q2) = shiftedQR(u, tol, iterCount);
    return {lam, q1 * q2};
}

} // namespace matoy
#endif