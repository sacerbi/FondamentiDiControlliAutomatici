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
// MAIN
// --------------------------
int main()
{
    // Dati
    double u1 = 1.0 / 8.0, u2 = 1.0 / 4.0, u3 = 1.0 / 2.0;
    double dx = 0.01;
    double x0 = -2, xf = 2;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i)
        x[i] = x0 + i * dx;
    // Calcolo il comportamento dei punti di equilibrio al variare di x
    std::vector<double> y1(N + 1);
    std::vector<double> y2(N + 1);
    std::vector<double> y3(N + 1);

    for (int i = 0; i <= N; i++)
    {
        y1[i] = x[i] * x[i] + x[i] + u1;
        y2[i] = x[i] * x[i] + x[i] + u2;
        y3[i] = x[i] * x[i] + x[i] + u3;
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
    fprintf(gp, "set title 'Variazione equilibrio al variare di ubar'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set xrange [-2:2]\n");
    fprintf(gp, "set yrange [-1:2]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'y(x), ubar = 0.125', '-' with lines title 'y(x), ubar = 0.25', '-' with lines title 'y(x), ubar = 0.5', 0 with lines dashtype 2 title 'y = 0'\n");

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
    // 3. Invio dei dati per la TERZA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], y3[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    return 0;
}
