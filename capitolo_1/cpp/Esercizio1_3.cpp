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
    double Mbar;
    double hbar;
    double kbar;
    double dFe;
    double sbar;
    double mu;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_anelloChiuso(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    (void)t; // unused
    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = -(p.hbar/p.Mbar) * x[1] - ((p.kbar + p.mu)/p.Mbar * x[0]) + ((p.kbar + p.mu)/p.Mbar) * p.sbar + p.dFe/p.Mbar;
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double Mbar = 1, hbar = 3, kbar = 1, sbar = 1, dFe = 0;
    std::vector<double> mus = {-kbar, -kbar-0.01};
    double s0 = 0.2, v0 = 0.2;
    double dt = 0.001;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Figura 1.3 - Transitorio posizione massa-molla sotto azione del controllore per valori di mu <= -kbar.
    // -----------------------------------------
    std::map<double, std::vector<double>> risultatiControllore;
    for (auto mu : mus) {
        Params p{Mbar, hbar, kbar, dFe, sbar, mu};
        double dummy = 0.0; // to match function signature
        auto traj = rungeKutta4(f_anelloChiuso, {s0, v0}, t0, tf, dt, p, dummy);
        std::vector<double> C;
        for (auto& x : traj) C.push_back(x[0]);
        risultatiControllore[mu] = C;
    }

    // Plot Figura 1.3
    plt::figure(1);
    for (auto& [mu, C] : risultatiControllore)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(mu)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della posizione - C.I. nominali");
    plt::legend();
    plt::grid(true);

    // -----------------------------------------
    // Figura 1.3b - Transitorio posizione massa-molla sotto azione del controllore per valori nominali e Fe ≠ Febar.
    // -----------------------------------------
    risultatiControllore.clear();
    dFe = 1.0;
    t0 = 0, tf = 100.0; // tempo finale più lungo per vedere l'effetto di dFe
    N = static_cast<int>((tf - t0) / dt);
    time.clear();
    time.resize(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    for (auto mu : mus) {
        Params p{Mbar, hbar, kbar, dFe, sbar, mu};
        double dummy = 0.0; // to match function signature
        auto traj = rungeKutta4(f_anelloChiuso, {s0, v0}, t0, tf, dt, p, dummy);
        std::vector<double> C;
        for (auto& x : traj) C.push_back(x[0]);
        risultatiControllore[mu] = C;
    }

    // Plot Figura 1.3b
    plt::figure(2);
    for (auto& [mu, C] : risultatiControllore)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(mu)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della posizione - C.I. nominali, Fe ≠ Febar");
    plt::legend();
    plt::grid(true);

    plt::show();

    return 0;
}
