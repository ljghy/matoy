#include <Matrix.hpp>
#include <cstdlib>
#include <iostream>
using namespace matoy;

using cp = std::complex<long double>;

template <typename Ty>
void printMat(const Matrix<Ty> &m, const char *comment = nullptr)
{
    if (comment)
        std::cout << comment << '\n';
    for (int i = 0; i < m.row(); ++i)
    {
        for (int j = 0; j < m.col(); ++j)
            std::cout << m(i, j) << ' ';
        std::cout << '\n';
    }
}

int main()
{
    Matrix<double> m{{1, 2, 3},
                     {4, 5, 6},
                     {7, 8, 8}};

    auto invm = rref(rowCat(m, eye<double>(3)))(Range::all, Range(3, 5));
    printMat(invm);
    auto p = qr(m);
    printMat(p[0], "Q =");
    printMat(p[1], "R =");
    printMat(p[0] * p[1]);

    Matrix<cp> cm{{1, cp(2, 1), 0, 3},
                  {4, 5, 0, 6},
                  {7, 8, 0, 8}};
    auto p2 = qr(cm);
    printMat(p2[0], "Q2 = ");
    printMat(p2[1], "R2 = ");
    printMat(p2[0] * p2[1] - cm);
}