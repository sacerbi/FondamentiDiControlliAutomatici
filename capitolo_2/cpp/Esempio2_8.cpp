#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <cstdio>
#include <thread>
#include <chrono>
#include "myDiffEquation.h"

// Gestione della compatibilità cross-platform per le pipe
#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

// --------------------------
// Parametri di sistema
// --------------------------
struct Params
{
    double J;
    double k;
};

double f_sign(double alpha)
{
    return alpha != 0 ? alpha / abs(alpha) : 0;
}

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_sistemaCentrifuga(double t, const std::vector<double> &x, const Params &p, double &u)
{
    (void)t; // unuser
    double dx = -(p.k / p.J) * x[0] * x[0] * f_sign(x[0]) + (1 / p.J) * u;
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main()
{
    // Dati
    double J = 1, k = 0.1, u = 0.5, w0 = 10;
    double dt = 0.01;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i)
        time[i] = t0 + i * dt;

    // -----------------------------------------
    // Andamento del sistema centrifuga
    // -----------------------------------------
    auto y = rungeKutta4(f_sistemaCentrifuga, {w0}, t0, tf, dt, Params{J, k}, u);
    std::vector<double> w;
    for (auto &v : y)
        w.push_back(v[0]);

    std::cout << "Simulation done. Plotting..." << std::endl;

    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return 1;
    }

    // ==========================================
    // FIGURA: Comportamento del sistema carrello-molla
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Comportamento del sistema centrifuga'\n");
    fprintf(gp, "set xlabel 'Tempo (t)'\n");
    fprintf(gp, "set ylabel 'Velocità angolare'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'w(t)'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < time.size(); ++i)
    {
        fprintf(gp, "%f %f\n", time[i], w[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    std::cout << "Starting second simulation." << std::endl;
    // Pausa di 1000 millisecondi per evitare la collisione sulla clipboard
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // Verifico con coppia motrice u(t) = 0
    u = 0;
    auto y1 = rungeKutta4(f_sistemaCentrifuga, {w0}, t0, tf, dt, Params{J, k}, u);
    std::vector<double> w1;
    for (auto &v : y1)
        w1.push_back(v[0]);
    std::vector<double> w1teorico(N + 1);
    for (int i = 0; i <= N; i++)
    {
        w1teorico[i] = (J * w0) / (J + k * w0 * f_sign(w0) * time[i]);
    }

    std::cout << "Second simulation done. Plotting..." << std::endl;

    // ==========================================
    // FIGURA: Comportamento del sistema centrifuga
    // ==========================================
    fprintf(gp, "set terminal qt 2 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Comportamento teorica del sistema centrifuga'\n");
    fprintf(gp, "set xlabel 'Tempo (t)'\n");
    fprintf(gp, "set ylabel 'Velocità angolare'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'w(t)', '-' with lines title 'w_t(t)'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < time.size(); ++i)
    {
        fprintf(gp, "%f %f\n", time[i], w1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < time.size(); ++i)
    {
        fprintf(gp, "%f %f\n", time[i], w1teorico[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    return 0;
}
