#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include "matplotlibcpp.h"
#include "myDiffEquation.h"

namespace plt = matplotlibcpp;

// --------------------------
// Parametri di sistema
// --------------------------
struct Params {
    double J;
    double k;
};

double f_sign(double alpha){
    return alpha != 0 ? alpha/abs(alpha) : 0;
}

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_sistemaCentrifuga(double t, const std::vector<double>& x, const Params& p, double& u) {
    (void)t; //unuser
    double dx = -(p.k/p.J)*x[0]*x[0]*f_sign(x[0]) + (1/p.J)*u;
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double J = 1, k = 0.1, u1 = 0.0, u2 = 0.5, w0 = 1;
    double dt = 0.01;
    double t0 = -5, tf = 5;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i) x[i] = t0 + i * dt;

    // Calcolo ybar all'equilibrio per entrambi gli input
    std::vector<double> w1(N+1);
    std::vector<double> w2(N+1);

    for (int i = 0; i<= N; i++){
        w1[i] = -(k/J) * x[i] * x[i] * f_sign(x[i]) + u1/J;
        w2[i] = -(k/J) * x[i] * x[i] * f_sign(x[i]) + u2/J;
    }
    // Plot Figura 2.4
    plt::figure(1);
    plt::plot(x, w1, {{"label", "w(t), u(t) = 0"}});
    plt::plot(x, w2, {{"label", "w(t), u(t) = 0.5"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Velocità angolare");
    plt::title("Comportamento del sistema centrifuga");
    plt::legend();
    plt::grid(true);
    
    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}
