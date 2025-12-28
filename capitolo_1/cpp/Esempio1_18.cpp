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
    double E1;
    double E2;
    double Fm1;
    double Fm2;
};

struct EditableParams {
    double Fm;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_anelloChiuso(double t, const std::vector<double>& x, const Params& p, EditableParams& ep) {
    (void)t; // unused
    std::vector<double> dx(2);
    double error = p.sbar - x[0];
    if (error >= p.E2) {
        ep.Fm = p.Fm2;
    }
    else if (error <= p.E1) {
        ep.Fm = p.Fm1;
    }
    else if (x[1] >= 0) {
        ep.Fm = p.Fm2;
    }
    else {
        ep.Fm = p.Fm1;
    }

    dx[0] = x[1];
    dx[1] = -(p.hbar/p.Mbar) * x[1] -(p.kbar/p.Mbar) * x[0] + (ep.Fm/p.Mbar) - p.Fe/p.Mbar;
    return dx;
}



// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double Mbar = 1, hbar = 3, kbar = 1, sbar = 1, Fe = 0;
    double E1 = -0.02, E2 = 0.02, Fm1 = -1, Fm2 = 2, Fm = 0;
    double dt = 0.001;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Figura 1.17 - Transitorio posizione massa-molla, C.I. nulle, sotto azione del controllore.
    // -----------------------------------------
    std::vector<double> transitorioPosizione;
    std::vector<double> transitorioFmotrice;
    Params p = Params{Mbar, hbar, kbar, Fe, sbar, E1, E2, Fm1, Fm2};
    EditableParams ep = EditableParams{Fm};

    auto traj = rungeKutta4(f_anelloChiuso, {0.0, 0.0}, t0, tf, dt, p, ep);
    // Estrae la posizione
    for (auto& x : traj) transitorioPosizione.push_back(x[0]);
    // Estrae la forza motrice
    for (int i = 0; i < traj.size(); ++i) {
        double error = sbar - traj[i][0];
        double Fm;
        if (error >= E2) {
            Fm = Fm2;
        }
        else if (error <= E1) {
            Fm = Fm1;
        }
        else if (traj[i][1] >= 0) {
            Fm = Fm2;
        }
        else {
            Fm = Fm1;
        }
        transitorioFmotrice.push_back(Fm);
    }

    // Plot Figura 1.17a - Transitorio posizione
    plt::figure(1);
    plt::plot(time, transitorioPosizione);
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della posizione - C.I. nulle");
    plt::legend();
    plt::grid(true);

    // Plot Figura 1.17b - Transitorio forza motrice
    plt::figure(2);
    plt::plot(time, transitorioFmotrice);
    plt::xlabel("Tempo (t)");
    plt::ylabel("Forza motrice Fm");
    plt::title("Transitorio della forza motrice - C.I. nulle");
    plt::legend();
    plt::grid(true);

    plt::show();

    return 0;
}
