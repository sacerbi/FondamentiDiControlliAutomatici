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
    double ubar; // placeholder se necessario
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_sistema(double t, const std::vector<double> &x, const Params &p, double &dummy)
{
    (void)t;                                 // unused
    double dx = x[0] * x[0] + x[0] + p.ubar; // ybar = x^2 + x + ubar
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main()
{
    // Dati
    double u1 = 0.125, u2 = 0.25;
    double dx = 0.01;
    double x0 = -2, xf = 2;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i)
        x[i] = x0 + i * dx;

    // Calcolo ybar all'equilibrio per entrambi gli input
    std::vector<double> y1(N + 1);
    std::vector<double> y2(N + 1);

    for (int i = 0; i <= N; i++)
    {
        y1[i] = x[i] * x[i] + x[i] + u1;
        y2[i] = x[i] * x[i] + x[i] + u2;
    }
    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return 1;
    }

    // ==========================================
    // FIGURA: Comportamento del sistema
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Variazione della funzione y(x) al variare di ubar'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y(x)'\n");
    fprintf(gp, "set xrange [%f:%f]\n", x0, xf);
    fprintf(gp, "set yrange [-1:2]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n");

    fprintf(gp, "plot '-' with lines title 'y(x), ubar = 0.125', '-' with lines title 'y(x), ubar = 0.25', 0 with lines title 'y=0'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], y1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], y2[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // Simulo per due condizioni di x(0) per ubar < 1/4
    // nella prima mi pongo prima del primo punto di stabilità
    // nella seconda mi pongo tra il primo e il secondo punto di stabilità
    double dt = 0.001;
    double t0 = 0.0, tf = 10.0;
    int M = static_cast<int>((tf - t0) / dt);
    std::vector<double> t(M + 1);
    for (int i = 0; i <= M; ++i)
        t[i] = t0 + i * dt;
    double x0_1 = -1.0; // prima del primo punto di stabilità
    double x0_2 = -0.3; // tra il primo e il secondo punto di stabilità
    double dummy = 0.0; // to match function signature
    auto y11 = rungeKutta4(f_sistema, {x0_1}, t0, tf, dt, Params{u1}, dummy);
    auto y12 = rungeKutta4(f_sistema, {x0_2}, t0, tf, dt, Params{u1}, dummy);

    std::vector<double> y11_vals;
    std::vector<double> y12_vals;
    for (auto &v : y11)
        y11_vals.push_back(v[0]);
    for (auto &v : y12)
        y12_vals.push_back(v[0]);

    // ==========================================
    // FIGURA: Evoluzione temporale
    // ==========================================
    fprintf(gp, "set terminal qt 2 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Evoluzione temporale di y(t) per due condizioni iniziali diverse'\n");
    fprintf(gp, "set xlabel 't'\n");
    fprintf(gp, "set ylabel 'y(t)'\n");
    fprintf(gp, "set xrange [%f:%f]\n", t0, tf);
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n");

    fprintf(gp, "plot '-' with lines title 'y(t), x0 = -1.0', '-' with lines title 'y(t), x0 = -0.3'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < t.size(); ++i)
    {
        fprintf(gp, "%f %f\n", t[i], y11_vals[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < t.size(); ++i)
    {
        fprintf(gp, "%f %f\n", t[i], y12_vals[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    return 0;
}
