#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

// --------------------------
// Parametri di sistema
// --------------------------
struct Params {
    double a;
    double b;
    double mu;
    double Cbar;
};

// --------------------------
// Equazione differenziale - Controllore proporzionale
// --------------------------
std::vector<double> f_prop(double t, const std::vector<double>& x, const Params& p) {
    (void)t; // unused
    std::vector<double> dx(1);
    dx[0] = -(p.a + p.b * p.mu) * x[0] + p.b * p.mu * p.Cbar;
    return dx;
}

// --------------------------
// Equazione differenziale - Controllore integrale
// --------------------------
std::vector<double> f_int(double t, const std::vector<double>& x, const Params& p) {
    (void)t;
    std::vector<double> dx(2);
    dx[0] = -p.a * x[0] + p.b * x[1];          // Cdot
    dx[1] = p.mu * p.Cbar - p.mu * x[0];       // Idot
    return dx;
}

// --------------------------
// Operatori utili
// --------------------------
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] + b[i];
    return r;
}
std::vector<double> operator*(const std::vector<double>& a, double s) {
    std::vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] * s;
    return r;
}

// --------------------------
// Runge–Kutta 4° ordine
// --------------------------
template <typename Func>
std::vector<std::vector<double>> rungeKutta4(Func f, std::vector<double> x0,
                                             double t0, double tf, double dt, const Params& p) {
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<std::vector<double>> result;
    result.reserve(N + 1);

    std::vector<double> x = x0;
    double t = t0;
    result.push_back(x);

    for (int i = 0; i < N; ++i) {
        auto k1 = f(t, x, p);
        auto k2 = f(t + dt / 2, x + k1 * (dt / 2), p);
        auto k3 = f(t + dt / 2, x + k2 * (dt / 2), p);
        auto k4 = f(t + dt, x + k3 * dt, p);

        x = x + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6);
        t += dt;
        result.push_back(x);
    }

    return result;
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double a = 1, b = 1, Cbar = 5;
    double dt = 0.01;

    // -----------------------------------------
    // Controllore proporzionale
    // -----------------------------------------
    std::vector<double> mus_prop = {1, 2, 4, 8, 16, 32};
    double t0 = 0, tf = 2.0;
    int N = static_cast<int>((tf - t0) / dt);

    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    std::map<double, std::vector<double>> Cprop;
    for (auto mu : mus_prop) {
        Params p{a, b, mu, Cbar};
        auto traj = rungeKutta4(f_prop, {0.0}, t0, tf, dt, p);
        std::vector<double> C;
        for (auto& x : traj) C.push_back(x[0]);
        Cprop[mu] = C;
    }

    // Plot proporzionale
    plt::figure(1);
    for (auto& [mu, C] : Cprop)
        plt::plot(time, C, {{"label", "\\u03BC = " + std::to_string(mu)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Stato C");
    plt::title("Evoluzione dello stato C - Controllore proporzionale");
    plt::legend();
    plt::grid(true);

    // -----------------------------------------
    // Controllore integrale
    // -----------------------------------------
    std::vector<double> mus_int = {0.1, 0.2, 0.4, 0.8, 1.6, 3.2};
    t0 = 0; tf = 12.0;
    N = static_cast<int>((tf - t0) / dt);

    std::vector<double> time2(N + 1);
    for (int i = 0; i <= N; ++i) time2[i] = t0 + i * dt;

    std::map<double, std::vector<double>> Cint;
    for (auto mu : mus_int) {
        Params p{a, b, mu, Cbar};
        auto traj = rungeKutta4(f_int, {0.0, 0.0}, t0, tf, dt, p);
        std::vector<double> C;
        for (auto& x : traj) C.push_back(x[0]);
        Cint[mu] = C;
    }

    // Plot integrale
    plt::figure(2);
    for (auto& [mu, C] : Cint)
        plt::plot(time2, C, {{"label", "\\u03BC = " + std::to_string(mu)}});
    plt::plot({0, 12}, {5, 5}, "k--");  // linea orizzontale Cbar
    plt::xlabel("Tempo (t)");
    plt::ylabel("Stato C");
    plt::title("Evoluzione dello stato C - Controllore integrale");
    plt::legend();
    plt::grid(true);

    plt::show();

    return 0;
}
