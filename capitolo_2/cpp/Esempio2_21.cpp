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
    double C;
    double L;
    double alpha;
    double beta;
};

struct Evoluzione {
    std::vector<double> x1;
    std::vector<double> x2;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_oscillatore(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    (void)t; // unused
    std::vector<double> dx(2);
    dx[0] = x[1]/p.C;
    dx[1] = - (3 * p.beta * x[0]*x[0] - p.alpha) * (x[1]/p.C) - (x[0]/p.L);
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double C = 1, L = 1, alpha = 1, beta = 1.0/3.0;
    std::vector<double> x1_0 = {-3, -2.5, -2, 3, 3, 2.5, 2, -0.01, 0.01};
    std::vector<double> x2_0 = {4, 4, 4, 4, -4, -4, -4, 0, 0};
    double dt = 0.001;
    double t0 = 0, tf = 40.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Andamento della tensione sul condensatore nel circuito RC e della tensione sul resistore
    // -----------------------------------------
    double dummy = 0.0; // to match function signature
    std::map<std::pair<double, double>, Evoluzione> evoluzioni;

    for (int i = 0; i < x1_0.size(); i++){
        std::vector<double> x0 = {x1_0[i], x2_0[i]};
        auto y = rungeKutta4(f_oscillatore, x0, t0, tf, dt, Params{C, L, alpha, beta}, dummy);
        std::vector<double> x1, x2;
        for (auto& v : y) {
            x1.push_back(v[0]);
            x2.push_back(v[1]);
        }
        evoluzioni[{x1_0[i], x2_0[i]}] = Evoluzione{x1, x2};
    }
    
    // Plot dei risultati
    plt::figure(1);
    for (auto& [chiave, evoluzione] : evoluzioni) {
        const auto& x1 = evoluzione.x1;
        const auto& x2 = evoluzione.x2;
        plt::plot(x1, x2, {{"label", "x1_0 = " + std::to_string(chiave.first) + ", x2_0 = " + std::to_string(chiave.second)}});
    }
    plt::xlabel("x1 (Tensione sul condensatore)");
    plt::ylabel("x2 (Corrente sul condensatore)");
    plt::title("Evoluzione dell'oscillatore di van der Pol");
    plt::legend();
    plt::grid(true);
    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}
