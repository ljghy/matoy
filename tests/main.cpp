#include <Matrix.hpp>
#include <cstdlib>
#include <iostream>
using namespace matoy;

using cp = std::complex<double>;

int main()
{
    Matrix<double> m{{-9, 0, 3},
                     {6, 5, 4},
                     {7, 8, -8}};

    std::cout << "inv(m) =\n"
              << inv(m);

    Matrix<double> q, r;
    std::tie(q, r) = qr(m);
    std::cout << "Q = \n"
              << q;
    std::cout << "R = \n"
              << r;

    Matrix<cp> cm{{1, cp(2, 1), 0, 3},
                  {4, 5, 0, 6},
                  {7, 8, 0, 8}};
    Matrix<cp> cq, cr;
    std::tie(cq, cr) = qr(cm);

    std::cout << "Q2 = \n"
              << cq;
    std::cout << "R2 = \n"
              << cr;

    Matrix<cp> m3{{cp(1, 1), cp(1, 2)}, {cp(-1, 0.5), cp(-2, 1)}};
    auto p3 = powit(m3);
    std::cout << p3 << '\n';
    std::cout << std::get<0>(eig(m3));

    Matrix<cp> ms{{1, 2, 3},
                  {2, 3, 4},
                  {3, 4, 1}};

    Matrix<cp> lam, ev;
    std::tie(lam, ev) = eig(ms);
    std::cout << real(lam)
              << real(ev);

    Matrix<cp> u, s, v;
    std::tie(u, s, v) = svd(rand<cp>(randnum<cp>, 3, 5));
    std::cout << "U = \n"
              << u;
    std::cout << "S = \n"
              << s;
    std::cout << "V = \n"
              << v;

    std::cout << "USVT = \n"
              << u * s * v.T();
}