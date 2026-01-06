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
    double dx = x[0] * x[0] + x[0] + p.ubar; // ybar = x^2 + x + ubar
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double u1 = 0.125, u2 = 0.25;
    double dx = 0.01;
    double x0 = -2, xf = 2;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i) x[i] = x0 + i * dx;

    // Calcolo ybar all'equilibrio per entrambi gli input
    std::vector<double> y1(N+1);
    std::vector<double> y2(N+1);

    for (int i = 0; i<= N; i++){
        y1[i] = x[i] * x[i] + x[i] + u1;
        y2[i] = x[i] * x[i] + x[i] + u2;
    }
    // Plot Figura 2.4
    plt::figure(1);
    plt::plot(x, y1, {{"label", "y(x), ubar = 0.125"}});
    plt::plot(x, y2, {{"label", "y(x), ubar = 0.25"}});
    plt::plot({x0, xf}, {0.0, 0.0}, {{"label", "y = 0"}, {"color", "black"}, {"linestyle", "-"}});
    plt::xlim(x0, xf);
    plt::ylim(-1, 2);
    plt::xlabel("x");
    plt::ylabel("y(x)");
    plt::title("Variazione della funzione y(x) al variare di ubar");
    plt::legend();
    plt::grid(true);

    // Simulo per due condizioni di x(0) per ubar < 1/4
    // nella prima mi pongo prima del primo punto di stabilità
    // nella seconda mi pongo tra il primo e il secondo punto di stabilità
    double dt = 0.001;
    double t0 = 0.0, tf = 10.0;
    int M = static_cast<int>((tf - t0) / dt);
    std::vector<double> t(M + 1);
    for (int i = 0; i <= M; ++i) t[i] = t0 + i * dt;
    double x0_1 = -1.0; // prima del primo punto di stabilità
    double x0_2 = -0.3;  // tra il primo e il secondo punto di stabilità
    double dummy = 0.0; // to match function signature
    auto y11 = rungeKutta4(f_sistema, {x0_1}, t0, tf, dt, Params{u1}, dummy);
    auto y12 = rungeKutta4(f_sistema, {x0_2}, t0, tf, dt, Params{u1}, dummy);

    std::vector<double> y11_vals;
    std::vector<double> y12_vals;
    for (auto& v : y11) y11_vals.push_back(v[0]);
    for (auto& v : y12) y12_vals.push_back(v[0]);

    plt::figure(2);
    plt::plot(t, y11_vals, {{"label", "y(t), x0 = -1.0"}});
    plt::plot(t, y12_vals, {{"label", "y(t), x0 = -0.3"}});
    plt::xlabel("t");
    plt::ylabel("y(t)");
    plt::title("Evoluzione temporale di y(t) per due condizioni iniziali diverse");
    plt::legend();
    plt::grid(true);
    
    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}
