#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include <cstdio>
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
    return alpha != 0 ? alpha / std::abs(alpha) : 0;
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
    double J = 1, k = 0.1, u1 = 0.0, u2 = 0.5, w0 = 1;
    double dt = 0.01;
    double t0 = -5, tf = 5;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i)
        x[i] = t0 + i * dt;

    // Calcolo ybar all'equilibrio per entrambi gli input
    std::vector<double> w1(N + 1);
    std::vector<double> w2(N + 1);

    for (int i = 0; i <= N; i++)
    {
        w1[i] = -(k / J) * x[i] * x[i] * f_sign(x[i]) + u1 / J;
        w2[i] = -(k / J) * x[i] * x[i] * f_sign(x[i]) + u2 / J;
    }

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

    fprintf(gp, "plot '-' with lines title 'w_1(t), u(t) = 0', '-' with lines title 'w_2(t), u(t) = 0.5'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], w1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], w2[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    return 0;
}
