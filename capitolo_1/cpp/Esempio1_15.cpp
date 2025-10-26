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
std::vector<double> f_anelloChiuso(double t, const std::vector<double>& x, const Params& p) {
    (void)t; // unused
    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = -(p.hbar/p.Mbar) * x[1] - ((p.kbar + p.mu)/p.Mbar * x[0]) + ((p.kbar + p.mu)/p.Mbar) * p.sbar + p.dFe/p.Mbar;
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
    double Mbar = 1, hbar = 3, kbar = 1, sbar = 1, dFe = 0;
    std::vector<double> mus_prop = {0.001, 0.5, 1.0, 5.0, 10.0, 20.0};
    double dt = 0.001;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Figura 1.12 - Transitorio posizione massa-molla, C.I. nulle, sotto azione del controllore.
    // -----------------------------------------
    std::map<double, std::vector<double>> risultatiControllore;
    for (auto mu : mus_prop) {
        Params p{Mbar, hbar, kbar, dFe, sbar, mu};
        auto traj = rungeKutta4(f_anelloChiuso, {0.0, 0.0}, t0, tf, dt, p);
        std::vector<double> C;
        for (auto& x : traj) C.push_back(x[0]);
        risultatiControllore[mu] = C;
    }

    // Plot Figura 1.12
    plt::figure(1);
    for (auto& [mu, C] : risultatiControllore)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(mu)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della posizione - C.I. nulle");
    plt::legend();
    plt::grid(true);

    // -----------------------------------------
    // Figura 1.13 - Transitorio posizione massa-molla sotto azione del controllore con sbar = 0.
    // a) C.I. non nulle
    // b) C.I. nulle e dFe != 0
    // -----------------------------------------
    risultatiControllore.clear();
    sbar = 0.0;
    double s0 = 1.0;
    double v0 = 1.0;
    for (auto mu : mus_prop) {
        Params p{Mbar, hbar, kbar, dFe, sbar, mu};
        auto traj = rungeKutta4(f_anelloChiuso, {s0, v0}, t0, tf, dt, p);
        std::vector<double> C;
        for (auto& x : traj) C.push_back(x[0]);
        risultatiControllore[mu] = C;
    }

    // Plot Figura 1.13a
    plt::figure(2);
    for (auto& [mu, C] : risultatiControllore)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(mu)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della posizione - C.I. non nulle");
    plt::legend();
    plt::grid(true);

    risultatiControllore.clear();
    dFe = 1.0;
    for (auto mu : mus_prop) {
        Params p{Mbar, hbar, kbar, dFe, sbar, mu};
        auto traj = rungeKutta4(f_anelloChiuso, {0.0, 0.0}, t0, tf, dt, p);
        std::vector<double> C;
        for (auto& x : traj) C.push_back(x[0]);
        risultatiControllore[mu] = C;
    }

    // Plot Figura 1.13b
    plt::figure(3);
    for (auto& [mu, C] : risultatiControllore)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(mu)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della posizione - C.I. nulle, Fe ≠ Febar");
    plt::legend();
    plt::grid(true);

    // -----------------------------------------
    // Figura 1.14 - Transitorio forza motricesotto azione del controllore con C.I. nulle.
    // -----------------------------------------
    std::map<double, std::vector<double>> risultatiFmotrice;
    sbar = 1.0;
    dFe = 0.0;
    for (auto mu : mus_prop) {
        Params p{Mbar, hbar, kbar, dFe, sbar, mu};
        auto traj = rungeKutta4(f_anelloChiuso, {0.0, 0.0}, t0, tf, dt, p);
        std::vector<double> C;
        for (auto& x : traj){
            double Fm = kbar * sbar + dFe + mu * (sbar - x[0]);
            C.push_back(Fm);
        }
        risultatiFmotrice[mu] = C;
    }

    // Plot Figura 1.14
    plt::figure(4);
    for (auto& [mu, C] : risultatiFmotrice)
        plt::plot(time, C, {{"label", "μ = " + std::to_string(mu)}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione s");
    plt::title("Transitorio della Fm - C.I. nulle");
    plt::legend();
    plt::grid(true);

    plt::show();

    return 0;
}
