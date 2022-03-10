#include <Matrix.hpp>
#include <cstdlib>
#include <iostream>
using namespace matoy;

using cp = std::complex<double>;

template <typename Ty>
void printMat(const Matrix<Ty> &m, const char *comment = nullptr)
{
    if (comment)
        std::cout << comment << '\n';
    for (int i = 0; i < m.row(); ++i)
    {
        for (int j = 0; j < m.col(); ++j)
            std::cout << m(i, j) << ' ';
        std::cout << "\n";
    }
}

int main()
{
    Matrix<double> m{{-9, 0, 3},
                     {6, 5, 4},
                     {7, 8, -8}};

    printMat(inv(m));

    Matrix<double> q, r;
    std::tie(q, r) = qr(m);
    printMat(q, "Q =");
    printMat(r, "R =");
    printMat(q * r);

    Matrix<cp> cm{{1, cp(2, 1), 0, 3},
                  {4, 5, 0, 6},
                  {7, 8, 0, 8}};
    Matrix<cp> cq, cr;
    std::tie(cq, cr) = qr(cm);
    printMat(cq, "Q2 = ");
    printMat(cr, "R2 = ");
    printMat(cq * cr);

    Matrix<cp> m3{{cp(1, 1), cp(1, 2)}, {cp(-1, 0.5), cp(-2, 1)}};
    auto p3 = powit(m3);
    std::cout << p3 << '\n';
    printMat(m3);
    printMat(std::get<0>(eig(m3)));

    Matrix<cp> ms{{1, 2, 3},
                  {2, 3, 4},
                  {3, 4, 1}};

    Matrix<cp> lam, ev;
    std::tie(lam, ev) = eig(ms);
    printMat(real(lam));
    printMat(real(ev));

    // Matrix<double> u, s, v;
    // std::tie(u, s, v) = svd(m({0, 1}, Range::all));
    // printMat(u, "U =");
    // printMat(s, "S =");
    // printMat(v, "V =");

}