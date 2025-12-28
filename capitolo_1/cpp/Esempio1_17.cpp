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
    double Fe;
    double sbar;
    double mu;
    double nu;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_anelloChiuso(double t, const std::vector<double>& x, const Params& p, double /*unused*/) {
    (void)t; // unused
    std::vector<double> dx(3);
    double error = p.sbar - x[0];
    dx[0] = x[1];
    dx[1] = -(p.hbar/p.Mbar) * x[1] -(p.kbar/p.Mbar) * x[0] + (p.mu/p.Mbar) * error + (p.nu/p.Mbar) * x[2] - p.Fe/p.Mbar;
    dx[2] = error;
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double Mbar = 1, hbar = 3, kbar = 1, sbar = 1, Fe = 0;
    std::vector<double> mus_prop = {1.0, 1.5, 3.0};
    std::vector<double> nus_prop = {0.5, 1.0};
    double dt = 0.001;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Figura 1.12 - Transitorio posizione massa-molla, C.I. nulle, sotto azione del controllore.
    // -----------------------------------------
    std::map<std::pair<double, double>, std::vector<double>> transitorioPosizione;
    std::map<std::pair<double, double>, std::vector<double>> transitorioFmotrice;

    for (auto mu : mus_prop) {
        for (auto nu : nus_prop) {
            Params p{Mbar, hbar, kbar, Fe, sbar, mu, nu};
            double dummy = 0.0; // parametro non usato nella funzione
            auto traj = rungeKutta4(f_anelloChiuso, {0.0, 0.0, 0.0}, t0, tf, dt, p, dummy);
            // Estrae la posizione
            std::vector<double> posizione;
            for (auto& x : traj) posizione.push_back(x[0]);
            transitorioPosizione[std::make_pair(mu, nu)] = posizione;
            // Estrae la forza motrice
            std::vector<double> Fmotrice;
            for (auto& x : traj) Fmotrice.push_back(mu * sbar - mu * x[0] + nu * x[2]);
            transitorioFmotrice[std::make_pair(mu, nu)] = Fmotrice;
        }
    }

    // Plot Figura 1.15a - Transitorio posizione
    plt::figure(1);
    for (auto& [params, C] : transitorioPosizione)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(params.first) + ", ν = " + std::to_string(params.second)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della posizione - C.I. nulle");
    plt::legend();
    plt::grid(true);

    // Plot Figura 1.15b - Transitorio forza motrice
    plt::figure(2);
    for (auto& [params, C] : transitorioFmotrice)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(params.first) + ", ν = " + std::to_string(params.second)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Forza motrice Fm");
    plt::title("Transitorio della forza motrice - C.I. nulle");
    plt::legend();
    plt::grid(true);

    plt::show();

    return 0;
}
