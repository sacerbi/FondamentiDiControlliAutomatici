// Per compilare con g++ ricordarsi di includere: -I./Eigen/ -std=c++20

#include <iostream>
#include <Dense>


int main() {
    // Definizione delle costanti (in Octave erano simboliche)
    double R = 1.0, C = 2.0;

    // Definizione delle matrici A, B, C
    // MatrixXd indica una matrice di dimensioni dinamiche di tipo double
    Eigen::MatrixXd Amat(2, 2);
    Amat << -1.0/(R*C), -1.0/(R*C),
            -1.0/(R*C), -1.0/(R*C);

    Eigen::MatrixXd Bmat(2, 1);
    Bmat << 1.0/(R*C),
            1.0/(R*C);

    Eigen::MatrixXd Cmat(1, 2);
    Cmat << 1, 0;

    // Matrice di raggiungibilità
    Eigen::MatrixXd Mr(2,2);
    Mr << Bmat, Amat*Bmat;
    std::cout << "Matrice di raggiungibilità:\n" << Mr << "\n\n";

    // Verifica di raggiungibilità
    Eigen::FullPivLU<Eigen::MatrixXd> luA(Mr);
    int rank = luA.rank();
    if (rank == Bmat.rows()){
        std::cout << "Completamente raggiungibile. Dimensione n:" << Bmat.rows() << ", rank: " << rank << "\n";
        return 0;
    }
    else
        std::cout << "Non completamente raggiungibile. Dimensione n:" << Bmat.rows() << ", rank: " << rank << "\n";

    // Se non completamente raggiungibile, scompongo
    // Matrice di trasformazione T
    Eigen::MatrixXd T1(2, 2);
    T1 << Mr(0,0), 0, Mr(1,0), 1;
    std::cout << "Initialized matrix T1\n" << T1 << "\n";
    auto T = T1.completeOrthogonalDecomposition().pseudoInverse();
    std::cout << "T matrix computed\n" << T << "\n"; 

    // Scomposizione
    Eigen::MatrixXd Ar = T * Amat * T1;
    Eigen::MatrixXd Br = T * Bmat;
    Eigen::MatrixXd Cr = Cmat * T1;

    double threshold = 0.0001;
    Ar = Ar.unaryExpr([threshold](double x) { 
        return (std::abs(x) < threshold) ? 0.0 : x; 
    });
    Br = Br.unaryExpr([threshold](double x) {
        return (std::abs(x) < threshold) ? 0.0 : x;
    });
    Cr = Cr.unaryExpr([threshold](double x) {
        return (std::abs(x) < threshold) ? 0.0 : x;
    });
    // Stampa dei risultati
    std::cout << "Matrice Ar:\n" << Ar << "\n\n";
    std::cout << "Matrice Br:\n" << Br << "\n\n";
    std::cout << "Matrice Cr:\n" << Cr << "\n";

    return 0;
}