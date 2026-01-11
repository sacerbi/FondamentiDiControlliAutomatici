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
    double ubar; // placeholder se necessario
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_sistema(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    (void)t; // unused
    double dx = abs(x[0]) + p.ubar;
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Considerando il sistema dinamico
    //   xdot(t) = |x(t)| + u(t)
    // Si calcoli il movimento dello stato per ubar=0 e x(0)=x0.
    // Poi per u(t) = ubar si determinino gli stati di equilibrio e le corrispondenti
    // proprietà di stabilità.
    double dt = 0.001;
    double t0 = 0.0, tf = 3.0;
    int M = static_cast<int>((tf - t0) / dt);
    std::vector<double> t(M + 1);
    for (int i = 0; i <= M; ++i) t[i] = t0 + i * dt;
    double ubar = 0.0, x0_1 = 1.0, x0_2 = -1.0, dummy = 0.0;
    auto dx1 = rungeKutta4(f_sistema, {x0_1}, t0, tf, dt, Params{ubar}, dummy);
    auto dx2 = rungeKutta4(f_sistema, {x0_2}, t0, tf, dt, Params{ubar}, dummy);
    std::vector<double> x1;
    std::vector<double> x2;
    for (auto& v : dx1) x1.push_back(v[0]);
    for (auto& v : dx2) x2.push_back(v[0]);
    plt::figure(1);
    plt::plot(t, x1, {{"label", "u = 0, x(0) = " + std::to_string(x0_1)}});
    plt::plot(t, x2, {{"label", "u = 0, x(0) = " + std::to_string(x0_2)}});
    plt::plot({t0, tf}, {0.0, 0.0}, {{"label", "x(t) = 0"}, {"color", "black"}, {"linestyle", "-"}});
    plt::xlabel("t");
    plt::ylabel("x(t)");
    plt::title("Esercizio 2.8 - Movimento stato");
    plt::legend();
    plt::grid(true);
    
    // Analisi degli stati
    double u1 = 1.0, u2 = 0.0, u3 = -1.0;
    double dx = 0.01;
    double x0 = -3, xf = 3;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i) x[i] = x0 + i * dx;

    // Calcolo ybar all'equilibrio per entrambi gli input
    std::vector<double> xdot1(N+1);
    std::vector<double> xdot2(N+1);
    std::vector<double> xdot3(N+1);

    for (int i = 0; i<= N; i++){
        xdot1[i] = abs(x[i]) + u1;
        xdot2[i] = abs(x[i]) + u2;
        xdot3[i] = abs(x[i]) + u3;
    }
    // Plot Figura 2.4
    plt::figure(2);
    plt::plot(x, xdot1, {{"label", "u = " + std::to_string(u1)}});
    plt::plot(x, xdot2, {{"label", "u = " + std::to_string(u2)}});
    plt::plot(x, xdot3, {{"label", "u = " + std::to_string(u3)}});
    plt::plot({x0, xf}, {0.0, 0.0}, {{"label", "y = 0"}, {"color", "black"}, {"linestyle", "-"}});
    plt::xlim(x0, xf);
    plt::xlabel("x");
    plt::ylabel("dx/dt");
    plt::title("Esercizio 2.8 - Studio equilibrio");
    plt::legend();
    plt::grid(true);
    
    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}
