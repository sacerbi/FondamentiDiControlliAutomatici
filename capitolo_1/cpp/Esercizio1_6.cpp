#include <iostream>
#include <vector>
#include <string>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double A = 0.5, hbar = 3, qbar = 2, Ebar = 0.5, h0 = 0, qprev = qbar;
    double dt = 0.001;
    double t0 = 0, tf = 3;
    // Vettori per i risultati
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;
    std::vector<double> h(N + 1);
    std::vector<double> q(N + 1);
    h[0] = h0; //Impongo la condizione iniziale

    // -----------------------------------------
    // Figura 1.6 - Evoluzione temporale del livello dell'acqua nel serbatoio
    // -----------------------------------------
    for (int i = 0; i < N; ++i) {
        double errore = hbar - h[i];
        if (errore >= Ebar)
            q[i] = qbar;
        else if (errore <= -Ebar)
            q[i] = -qbar;
        else
            q[i] = qprev;
        qprev = q[i];
        // Calcolo della derivata di h(t)
        double hdot = q[i] / A;
        // Calcolo del nuovo valore di h(t+dt)
        h[i + 1] = h[i] + hdot * dt;
    }

    
    // Plot Figura 1.6 - Livello dell'acqua nel serbatoio
    plt::figure(1);
    plt::plot(time, h, {{"label", "h(t)"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Livello dell'acqua h");
    plt::title("Evoluzione temporale del livello dell'acqua nel serbatoio");
    plt::legend();
    plt::grid(true);

    // Plot Figura 1.6b - Portata
    plt::figure(2);
    plt::plot(time, q, {{"label", "q(t)"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Portata q");
    plt::title("Evoluzione temporale della portata nel serbatoio");
    plt::legend();
    plt::grid(true);

    plt::show();

    return 0;
}
