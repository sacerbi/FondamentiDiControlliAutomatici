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
    double dx = x[0] > 0 ? (p.ubar - 2*(1-exp(-x[0]))) : (p.ubar + x[0]);
    return {dx};
}

template<typename Func>
std::pair<std::vector<double>,std::vector<double>> evolviSistema(Func sistema, double ti, double tf, double dt, double x0, const Params& p){
    int M = static_cast<int>((tf - ti) / dt);
    std::vector<double> t(M + 1);
    for (int i = 0; i <= M; ++i) t[i] = ti + i * dt;
    double dummy = 0.0; // to match function signature
    auto y = rungeKutta4(sistema, {x0}, ti, tf, dt, p, dummy);
    std::vector<double> y_vals;
    for (auto& v : y) y_vals.push_back(v[0]);
    return {t, y_vals};
}

void esercizio2_7_1(std::vector<double> x) {
    // Dati
    double u1 = -1.0, u2 = 0.0, u3 = 1.0;
    std::vector<double> xp1(x.size());
    std::vector<double> xp2(x.size());
    std::vector<double> xp3(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        xp1[i] = u1 * (x[i]*x[i]*x[i]);
        xp2[i] = u2 * (x[i]*x[i]*x[i]);
        xp3[i] = u3 * (x[i]*x[i]*x[i]);
    }
    // Plot Figura 2.7.1
    plt::figure(1);
    plt::plot(x, xp1, {{"label", "ubar = -1.0"}});
    plt::plot(x, xp2, {{"label", "ubar = 0.0"}});
    plt::plot(x, xp3, {{"label", "ubar = 1.0"}});
    plt::xlim(x.front(), x.back());
    plt::ylim(-1, 2);
    plt::xlabel("x");
    plt::ylabel("dx/dt");
    plt::title("Esercizio 2.7.1");
    plt::legend();
    plt::grid(true);

    plt::show();
}

void esercizio2_7_2(std::vector<double> x) {
    // Dati
    double u1 = -1.0, u2 = 0.0, u3 = 1.0, u4 = 2.1;
    std::vector<double> xp1(x.size());
    std::vector<double> xp2(x.size());
    std::vector<double> xp3(x.size());
    std::vector<double> xp4(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        if(x[i] <= 0){
            xp1[i] = u1 + x[i];
            xp2[i] = u2 + x[i];
            xp3[i] = u3 + x[i];
            xp4[i] = u4 + x[i];
        }
        else{
            xp1[i] = u1 - 2*(1-exp(-x[i]));
            xp2[i] = u2 - 2*(1-exp(-x[i]));
            xp3[i] = u3 - 2*(1-exp(-x[i]));
            xp4[i] = u4 - 2*(1-exp(-x[i]));
        }
    }
    // Plot Figura 2.7.2
    plt::figure(2);
    plt::plot(x, xp1, {{"label", "ubar = -1.0"}});
    plt::plot(x, xp2, {{"label", "ubar = 0.0"}});
    plt::plot(x, xp3, {{"label", "ubar = 1.0"}});
    plt::plot(x, xp4, {{"label", "ubar = 2.1"}});
    plt::plot({x.front(), x.back()}, {0.0, 0.0}, {{"label", "y = 0"}, {"color", "black"}, {"linestyle", "-"}});
    plt::xlim(x.front(), x.back());
    plt::xlabel("x");
    plt::ylabel("dx/dt");
    plt::title("Esercizio 2.7.2 - Equilibrio");
    plt::legend();
    plt::grid(true);

    plt::show();

    auto [t1, xdot1] = evolviSistema(f_sistema, 0.0, 10.0, 0.001, -0.9, Params{1.0});
    auto [t2, xdot2] = evolviSistema(f_sistema, 0.0, 10.0, 0.001, 3.0, Params{1.0});
    auto xbar = -log(0.5);

    plt::figure(3);
    plt::plot(t1, xdot1, {{"label", "u(t) = 1, x(0) = -0.9"}});
    plt::plot(t1, xdot2, {{"label", "u(t) = 1, x(0) = 3.0"}});
    plt::plot({t1.front(), t1.back()}, {xbar, xbar}, {{"label", "x = " + std::to_string(xbar)}, {"color", "black"}, {"linestyle", "-"}});
    plt::xlabel("t");
    plt::ylabel("x(t)");
    plt::title("Esercizio 2.7.2 - Stabilità");
    plt::legend();
    plt::grid(true);

    plt::show();
}

void esercizio2_7_3(std::vector<double> x){
    std::map<std::pair<double,double>, std::vector<double>> risultati;
    // Caso ubar = 0, diversi valori di alpha
    double ubar = 0;
    std::vector<double> alphas = {-1.0, 0.0, 1.0};
    for(auto alpha : alphas){
        std::vector<double> xdot(x.size());
        for(int i = 0; i<xdot.size(); i++){
            xdot[i] = alpha * x[i] - ubar * ubar * x[i] * x[i] * x[i];
        }
        risultati[{ubar, alpha}] = xdot;
    }
    // Caso ubar != 0, alpha diversi
    ubar = 1.0;
    for(auto alpha : alphas){
        std::vector<double> xdot(x.size());
        for(int i = 0; i<xdot.size(); i++){
            xdot[i] = alpha * x[i] - ubar * ubar * x[i] * x[i] * x[i];
        }
        risultati[{ubar, alpha}] = xdot;
    }

    // Plot dei risultati
    plt::figure(4);
    for (auto& [chiave, evoluzione] : risultati) {
        plt::plot(x, evoluzione, {{"label", "u = " + std::to_string(chiave.first) + ", α = " + std::to_string(chiave.second)}});
    }
    plt::plot({x.front(), x.back()}, {0.0, 0.0}, {{"label", "y = 0"}, {"color", "black"}, {"linestyle", "-"}});
    plt::xlabel("x");
    plt::ylabel("dx/dt");
    plt::title("Esercizio 2.7.3");
    plt::legend();
    plt::grid(true);
    plt::show();
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double u1 = 0.125, u2 = 0.25;
    double dx = 0.01;
    double x0 = -3, xf = 3;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i) x[i] = x0 + i * dx;

    // Punto 1
    esercizio2_7_1(x);

    // Punto 2
    esercizio2_7_2(x);

    // Punto 3
    esercizio2_7_3(x);

    plt::detail::_interpreter::kill();

    return 0;
}
