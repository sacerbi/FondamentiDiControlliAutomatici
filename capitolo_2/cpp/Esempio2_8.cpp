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
    double J;
    double k;
};

double f_sign(double alpha){
    return alpha != 0 ? alpha/abs(alpha) : 0;
}

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_sistemaCentrifuga(double t, const std::vector<double>& x, const Params& p, double& u) {
    (void)t; //unuser
    double dx = -(p.k/p.J)*x[0]*x[0]*f_sign(x[0]) + (1/p.J)*u;
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double J = 1, k = 0.1, u = 0.5, w0 = 10;
    double dt = 0.01;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Andamento della tensione sul condensatore nel circuito RC e della tensione sul resistore
    // -----------------------------------------
    auto y = rungeKutta4(f_sistemaCentrifuga, {w0}, t0, tf, dt, Params{J, k}, u);
    std::vector<double> w;
    for(auto& v : y) w.push_back(v[0]);
    // Plot Figura 2.4
    plt::figure(1);
    plt::plot(time, w, {{"label", "w(t)"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Velocità angolare");
    plt::title("Comportamento del sistema centrifuga");
    plt::legend();
    plt::grid(true);

    // Verifico con coppia motrice u(t) = 0
    u = 0;
    auto y1 = rungeKutta4(f_sistemaCentrifuga, {w0}, t0, tf, dt, Params{J, k}, u);
    std::vector<double> w1;
    for(auto& v : y1) w1.push_back(v[0]);
    std::vector<double> w1teorico(N+1);
    for (int i = 0; i <= N; i++){
        w1teorico[i] = (J*w0)/(J + k*w0*f_sign(w0)*time[i]);
    }
    plt::figure(2);
    plt::plot(time, w1, {{"label", "w(t)"}});
    plt::plot(time, w1teorico, {{"label", "w(t) teorico"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Velocità angolare");
    plt::title("Comportamento del sistema centrifuga con coppia motrice nulla");
    plt::legend();
    plt::grid(true);
    
    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}
