#include <iostream>
#include "Eigen/Dense"

int main()
{
    // Definizione delle costanti (in Octave erano simboliche)
    double k = 1.0, M1 = 2.0, M2 = 2.0;

    // Definizione delle matrici A, B, C
    // MatrixXd indica una matrice di dimensioni dinamiche di tipo double
    Eigen::MatrixXd A(4, 4);
    A << 0, 0, 1, 0,
        0, 0, 0, 1,
        -k / M1, k / M1, 0, 0,
        k / M2, -k / M2, 0, 0;

    Eigen::MatrixXd B(4, 2);
    B << 0, 0,
        0, 0,
        1 / M1, 0,
        0, 1 / M2;

    Eigen::MatrixXd C(1, 4);
    C << 1, 0, 0, 0;

    // Matrice di trasformazione T
    Eigen::MatrixXd T(4, 4);
    T << 1, 0, 0, 0,
        -1, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, -1, 1;

    // Calcolo delle matrici trasformate
    // pinv() in Octave è la pseudo-inversa.
    // In Eigen si usa completeOrthogonalDecomposition per la pseudo-inversa
    Eigen::MatrixXd T_pinv = T.completeOrthogonalDecomposition().pseudoInverse();

    Eigen::MatrixXd At = T * A * T_pinv;
    Eigen::MatrixXd Bt = T * B;
    Eigen::MatrixXd Ct = C * T_pinv;

    double threshold = 0.0001;
    At = At.unaryExpr([threshold](double x)
                      { return (std::abs(x) < threshold) ? 0.0 : x; });
    Bt = Bt.unaryExpr([threshold](double x)
                      { return (std::abs(x) < threshold) ? 0.0 : x; });
    Ct = Ct.unaryExpr([threshold](double x)
                      { return (std::abs(x) < threshold) ? 0.0 : x; });
    // Stampa dei risultati
    std::cout << "Matrice At:\n"
              << At << "\n\n";
    std::cout << "Matrice Bt:\n"
              << Bt << "\n\n";
    std::cout << "Matrice Ct:\n"
              << Ct << "\n";

    return 0;
}