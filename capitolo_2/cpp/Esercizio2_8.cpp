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
    (void)t; // unused
    double dx = std::abs(x[0]) + p.ubar;
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main()
{
    // Considerando il sistema dinamico
    //   xdot(t) = |x(t)| + u(t)
    // Si calcoli il movimento dello stato per ubar=0 e x(0)=x0.
    // Poi per u(t) = ubar si determinino gli stati di equilibrio e le corrispondenti
    // proprietà di stabilità.
    double dt = 0.001;
    double t0 = 0.0, tf = 3.0;
    int M = static_cast<int>((tf - t0) / dt);
    std::vector<double> t(M + 1);
    for (int i = 0; i <= M; ++i)
        t[i] = t0 + i * dt;
    double ubar = 0.0, x0_1 = 1.0, x0_2 = -1.0, dummy = 0.0;
    auto dx1 = rungeKutta4(f_sistema, {x0_1}, t0, tf, dt, Params{ubar}, dummy);
    auto dx2 = rungeKutta4(f_sistema, {x0_2}, t0, tf, dt, Params{ubar}, dummy);
    std::vector<double> x1;
    std::vector<double> x2;
    for (auto &v : dx1)
        x1.push_back(v[0]);
    for (auto &v : dx2)
        x2.push_back(v[0]);
    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return 1;
    }

    // ==========================================
    // FIGURA: Esercizio 2.8 - Movimento stato
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Esercizio 2.8 - Movimento stato'\n");
    fprintf(gp, "set xlabel 't'\n");
    fprintf(gp, "set ylabel 'x(t)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'u = 0, x(0) = %f', '-' with lines title 'u = 0, x(0) = %f', 0 with lines title 'x(t) = 0'\n", x0_1, x0_2);

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < t.size(); ++i)
    {
        fprintf(gp, "%f %f\n", t[i], x1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < t.size(); ++i)
    {
        fprintf(gp, "%f %f\n", t[i], x2[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    std::cout << "Starting second simulation." << std::endl;
    // Pausa di 1000 millisecondi per evitare la collisione sulla clipboard
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    // Analisi degli stati
    double u1 = 1.0, u2 = 0.0, u3 = -1.0;
    double dx = 0.01;
    double x0 = -3, xf = 3;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i)
        x[i] = x0 + i * dx;

    // Calcolo ybar all'equilibrio per entrambi gli input
    std::vector<double> xdot1(N + 1);
    std::vector<double> xdot2(N + 1);
    std::vector<double> xdot3(N + 1);

    for (int i = 0; i <= N; i++)
    {
        xdot1[i] = std::abs(x[i]) + u1;
        xdot2[i] = std::abs(x[i]) + u2;
        xdot3[i] = std::abs(x[i]) + u3;
    }
    // ==========================================
    // FIGURA: Esercizio 2.8 - Studio equilibrio
    // ==========================================
    fprintf(gp, "set terminal qt 2 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Esercizio 2.8 - Studio equilibrio'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'dx / dt'\n");
    fprintf(gp, "set xrange [%f:%f]\n", x0, xf);
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'u = %f', '-' with lines title 'u = %f', '-' with lines title 'u = %f', 0 with lines title 'x(t) = 0'\n", u1, u2, u3);

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xdot1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xdot2[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
                        // 3. Invio dei dati per la TERZA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xdot3[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    return 0;
}
