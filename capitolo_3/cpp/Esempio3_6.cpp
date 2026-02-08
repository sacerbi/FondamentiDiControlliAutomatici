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
    double k, h, M;
    std::vector<double> u;
    double dt;
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_massSpringDamper(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    int idx = static_cast<int>(t / p.dt);
    double ut = idx < p.u.size() ? p.u[idx] : p.u.back();
    
    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = -(p.k/p.M)*x[0] - (p.h/p.M)*x[1] + (1/p.M)*ut; // Includo nel calcolo la Fext
    return dx;
}

std::vector<double> f_caso1(const std::vector<double>& t, const std::vector<double>& x0, const Params& p){
    double disc = (p.h * p.h) / (4 * p.M * p.M) - p.k / p.M;
    double s1 = -p.h / (2 * p.M) + std::sqrt(disc);
    double s2 = -p.h / (2 * p.M) - std::sqrt(disc);
    std::vector<double> y_teorica;
    for (int i = 0; i < t.size(); i++) {
        double yl = (1.0/(s2-s1)) * ((s2*x0[0]-x0[1])*exp(s1*t[i]) - (s1*x0[0]-x0[1])*exp(s2*t[i]));
        double yf = (1.0/(p.M*(s1-s2))) * ((exp(s1*t[i])-1)/s1 - (exp(s2*t[i])-1)/s2);
        y_teorica.push_back(yl + yf);
    }
    return y_teorica;
}

std::vector<double> f_caso2(const std::vector<double>& t, const std::vector<double>& x0, const Params& p){
    double sigma = -p.h/(2*p.M);
    double omega = std::sqrt(p.k/p.M - (p.h*p.h)/(4*p.M*p.M));
    double delta = std::sqrt(sigma*sigma + omega*omega);
    double gamma = std::acos(sigma/delta);
    std::vector<double> y_teorica;
    for (int i = 0; i < t.size(); i++) {
        double yl = exp(sigma*t[i]) * (-(delta/omega)*std::sin(omega*t[i]-gamma)*x0[0] + (1/omega)*std::sin(omega*t[i])*x0[1]);
        double yf = 1/p.M * (1.0/(delta*delta) * (1+delta/omega * exp(sigma*t[i]) * std::sin(omega*t[i]-gamma)));
        y_teorica.push_back(yl+yf);
    }
    return y_teorica;
}

std::vector<double> f_caso3(const std::vector<double>& t, const std::vector<double>& x0, const Params& p){
    double s0 = -p.h/(2*p.M);
    std::vector<double> y_teorica;
    for (int i = 0; i < t.size(); i++) {
        double yl = exp(s0*t[i]) * (x0[0] - (s0*x0[0] - x0[1])*t[i]);
        double yf = 1.0/p.M * (1.0/(s0*s0) - exp(s0*t[i])/(s0*s0) + (t[i]*exp(s0*t[i]))/s0);
        y_teorica.push_back(yl+yf);
    }
    return y_teorica;
}

double rms(const std::vector<double>& s1, const std::vector<double>& s2, double n){
    for (int i = 0; i < s1.size(); ++i) {
        n += std::pow(s1[i] - s2[i], 2);
    }
    n = std::sqrt(n / (s1.size() + 1));
    return n;
}

int main() {
    // Configurazione Caso 1: h^2 > 4Mk
    double dt = 0.001;
    double ti = 0.0, tf = 50.0;
    std::vector<double> x0 = {2.0, 1.0};
    int N = static_cast<int>((tf-ti)/dt);
    std::vector<double> time(N+1);
    for (int i = 0; i <= N; ++i) time[i] = ti + i*dt;
    std::vector<double> u(N+1);
    for (int i = 0; i <= N; ++i) u[i] = 1.0;
    Params p1 = {1.0, 3.0, 2.0, u, dt};
    
    // Simulazione del sistema
    double dummy = 0.0;
    auto dx1 = rungeKutta4(f_massSpringDamper, x0, ti, tf, dt, p1, dummy);

    std::vector<double> y1;
    for (auto& x : dx1) y1.push_back(x[0]);

    // Confronto del sistema con la soluzione teorica
    std::vector<double> y1_teorica = f_caso1(time, x0, p1);

    // Plot Caso 1
    plt::figure(1);
    plt::plot(time, y1, {{"label", "y"}});
    plt::plot(time, y1_teorica, {{"label", "y_t"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione (m)");
    plt::title("Comportamento del sistema k,h,M");
    plt::legend();
    plt::grid(true);
    // Calcolo dell'errore
    double error1 = rms(y1, y1_teorica, 0.0);
    std::cout << "Errore RMS y: " << error1 << std::endl;

    // Configurazione Caso 2: h^2 < 4Mk
    Params p2 = {1.0, 1.0, 2.0, u, dt};
    auto dx2 = rungeKutta4(f_massSpringDamper, x0, ti, tf, dt, p2, dummy);

    std::vector<double> y2;
    for (auto& x : dx2) y2.push_back(x[0]);

    // Confronto del sistema con la soluzione teorica
    std::vector<double> y2_teorica = f_caso2(time, x0, p2);

    // Plot Caso 2
    plt::figure(2);
    plt::plot(time, y2, {{"label", "y"}});
    plt::plot(time, y2_teorica, {{"label", "y_t"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione (m)");
    plt::title("Comportamento del sistema k,h,M");
    plt::legend();
    plt::grid(true);
    // Calcolo dell'errore
    double error2 = rms(y2, y2_teorica, 0.0);
    std::cout << "Errore RMS y: " << error2 << std::endl;


    // Configurazione Caso 1: h^2 == 4Mk
    Params p3 = {2.0, 4.0, 2.0, u, dt};
    auto dx3 = rungeKutta4(f_massSpringDamper, x0, ti, tf, dt, p3, dummy);

    std::vector<double> y3;
    for (auto& x : dx3) y3.push_back(x[0]);

    // Confronto del sistema con la soluzione teorica
    std::vector<double> y3_teorica = f_caso3(time, x0, p3);

    // Plot Caso 3
    plt::figure(3);
    plt::plot(time, y3, {{"label", "y"}});
    plt::plot(time, y3_teorica, {{"label", "y_t"}});
    plt::xlabel("Tempo (t)");
    plt::ylabel("Posizione (m)");
    plt::title("Comportamento del sistema k,h,M");
    plt::legend();
    plt::grid(true);
    // Calcolo dell'errore
    double error3 = rms(y3, y3_teorica, 0.0);
    std::cout << "Errore RMS y: " << error3 << std::endl;

    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}