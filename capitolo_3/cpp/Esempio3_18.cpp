#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <cmath>
#include "myDiffEquation.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// --------------------------
// Parametri del sistema
// --------------------------
struct Params{
    double g, k, L, M;
    std::vector<double> u;
    double dt;
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_pendolo(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    int idx = static_cast<int>(t / p.dt);
    double ut = idx < p.u.size() ? p.u[idx] : p.u.back();
    
    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = -(p.g/p.L)*sin(x[0]) - (p.k/(p.M*p.L*p.L))*x[1] + 1/(p.M*p.L*p.L) * ut;
    return dx;
}

int main() {
    // Dati
    double M = 10.0, k = 2.0, g=9.81, L=0.5;
    double dt = 0.001;
    double ti = 0.0, tf = 10.0;
    std::vector<double> x0 = {M_PI / 2, 0.0};
    int N = static_cast<int>((tf-ti)/dt);
    std::vector<double> time(N+1);
    for (int i = 0; i <= N; ++i) time[i] = ti + i*dt;
    std::vector<double> u(N+1);
    for (int i = 0; i <= N; ++i) u[i] = 0.0;
    Params p1 = {g, k, L, M, u, dt};
    
    // Simulazione del sistema
    double dummy = 0.0;
    auto dx1 = rungeKutta4(f_pendolo, x0, ti, tf, dt, p1, dummy);

    std::vector<double> x1, x2, y;
    for (auto& x : dx1) {
        x1.push_back(x[0]);
        x2.push_back(x[1]);
        double x2_val = (1/2)*M*L*L*x[1]*x[1];
        double x1_val = -M*g*L*cos(x[0]);
        y.push_back(x2_val+x1_val);
    }

    // Plot Caso 1
    plt::figure(1);
    plt::plot(time, y, {{"label", "E"}});
    plt::plot(time, x1, {{"label", "pos"}});
    plt::plot(time, x2, {{"label", "vel"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione/Velocità/Energia");
    plt::title("Comportamento del sistema pendolo");
    plt::legend();
    plt::grid(true);

    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}