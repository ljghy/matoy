#ifndef MATOY_IMPL_VECTOR_INL_
#define MATOY_IMPL_VECTOR_INL_

namespace matoy
{

template <typename Ty>
inline Ty sqrNorm(const Matrix<Ty> &vec)
{
    assert(vec.col() == 1);
    return vec.T() * vec;
}

template <typename Ty>
inline Ty sqrNorm(const Matrix<std::complex<Ty>> &vec)
{
    assert(vec.col() == 1);
    Ty s = 0.0;
    for (index_t i = 0; i < vec.row(); ++i)
        s += std::real(std::conj(vec(i, 0)) * vec(i, 0));
    return s;
}

template <typename Ty>
inline Ty norm(const Matrix<Ty> &vec)
{
    assert(vec.col() == 1);
    return std::sqrt(vec.T() * vec);
}

template <typename Ty>
inline Ty norm(const Matrix<std::complex<Ty>> &vec)
{
    assert(vec.col() == 1);
    return std::sqrt(sqrNorm(vec));
}

template <typename Ty>
inline Matrix<Ty> orthoBasis(index_t dim, index_t axis, const Ty &len = 1.0)
{
    Matrix<Ty> ret(dim, 1);
    for (index_t i = 0; i < dim; ++i)
        ret(i, 0) = 0.0;
    ret(axis, 0) = len;
    return ret;
}

template <typename Ty>
inline Matrix<Ty> normalized(const Matrix<Ty> &vec)
{
    return vec / norm(vec);
}

template <typename Ty>
inline void normalize(Matrix<Ty> &vec)
{
    vec /= norm(vec);
}

template <typename Ty>
inline Matrix<Ty> diag(const Matrix<Ty> &mat)
{
    index_t n = std::min(mat.row(), mat.col());
    Matrix<Ty> ret(n, 1);
    for (index_t i = 0; i < n; ++i)
        ret(i, 0) = mat(i, i);
    return ret;
}

} // namespace matoy

#endif