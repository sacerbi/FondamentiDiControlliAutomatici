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
    double M;
    double Fm;
    double kt0;
    double alpha;
    double t0;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_sistemaCarrelloMolla(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    double u = p.Fm;
    double Fe = p.kt0 * exp(-p.alpha*(t-p.t0)) * x[0];
    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = (u - Fe)/p.M;
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double Fm = 5, kt0 = 10, alpha = 2, s0 = 0, s0dot = 0, M = 10;
    double dt = 0.01;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Andamento della tensione sul condensatore nel circuito RC e della tensione sul resistore
    // -----------------------------------------
    double dummy = 0.0; // to match function signature
    auto y = rungeKutta4(f_sistemaCarrelloMolla, {s0, s0dot}, t0, tf, dt, Params{M, Fm, kt0, alpha, t0}, dummy);

    std::vector<double> s, sdot, et;
    for (auto& v : y) {
        s.push_back(v[0]);
        sdot.push_back(v[1]);
    }
    for (int i = 0; i <= N; ++i){
        double Ue = (kt0 * exp(-alpha*(time[i]-t0)) * s[i]*s[i])/2;
        double Up = (M*sdot[i]*sdot[i])/2;
        et.push_back(Ue+Up);
    }

    // Plot Figura 2.4
    plt::figure(1);
    plt::plot(time, s, {{"label", "s(t)"}});
    plt::plot(time, sdot, {{"label", "sdot(t)"}});
    plt::plot(time, et, {{"label", "ET(t)"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Spostamento / Velocità / Energia Totale");
    plt::title("Comportamento del sistema carrello-molla");
    plt::legend();
    plt::grid(true);
    
    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}
