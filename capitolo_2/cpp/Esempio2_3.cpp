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
    double R;
    double C;
    double Vg;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_circuitoRC(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    (void)t; // unused
    double dx = (1/(p.R * p.C)) * (p.Vg - x[0]);
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double R = 1, C = 2, Vg = 5, x0 = 1;
    double dt = 0.001;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Andamento della tensione sul condensatore nel circuito RC e della tensione sul resistore
    // -----------------------------------------
    double dummy = 0.0; // to match function signature
    auto y = rungeKutta4(f_circuitoRC, {x0}, t0, tf, dt, Params{R, C, Vg}, dummy);

    std::vector<double> Vc, Vr;
    for (auto& v : y) Vc.push_back(v[0]);
    for (auto& v : y) Vr.push_back(Vg - v[0]);

    // Plot Figura 1.12
    plt::figure(1);
    plt::plot(time, Vc, {{"label", "Vc"}});
    plt::plot(time, Vr, {{"label", "Vr"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Tensione (V)");
    plt::title("Comportamento del circuito RC");
    plt::legend();
    plt::grid(true);

    plt::show();

    return 0;
}
