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
    double C;
    double L;
    double alpha;
    double beta;
};

struct Evoluzione
{
    std::vector<double> x1;
    std::vector<double> x2;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_oscillatore(double t, const std::vector<double> &x, const Params &p, double &dummy)
{
    (void)t; // unused
    std::vector<double> dx(2);
    dx[0] = x[1] / p.C;
    dx[1] = -(3 * p.beta * x[0] * x[0] - p.alpha) * (x[1] / p.C) - (x[0] / p.L);
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main()
{
    // Dati
    double C = 1, L = 1, alpha = 1, beta = 1.0 / 3.0;
    std::vector<double> x1_0 = {-3, -2.5, -2, 3, 3, 2.5, 2, -0.01, 0.01};
    std::vector<double> x2_0 = {4, 4, 4, 4, -4, -4, -4, 0, 0};
    double dt = 0.01;
    double t0 = 0, tf = 40.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i)
        time[i] = t0 + i * dt;

    // -----------------------------------------
    // Andamento della tensione sul condensatore nel circuito RC e della tensione sul resistore
    // -----------------------------------------
    double dummy = 0.0; // to match function signature
    std::map<std::pair<double, double>, Evoluzione> evoluzioni;

    for (int i = 0; i < x1_0.size(); i++)
    {
        std::vector<double> x0 = {x1_0[i], x2_0[i]};
        auto y = rungeKutta4(f_oscillatore, x0, t0, tf, dt, Params{C, L, alpha, beta}, dummy);
        std::vector<double> x1, x2;
        for (auto &v : y)
        {
            x1.push_back(v[0]);
            x2.push_back(v[1]);
        }
        evoluzioni[{x1_0[i], x2_0[i]}] = Evoluzione{x1, x2};
    }

    // Plot dei risultati
    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return 1;
    }

    // ==========================================
    // FIGURA: Evoluzione dell'oscillatore di van der Pol
    // ==========================================
    // fprintf(gp, "set terminal qt 1 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title \"Evoluzione dell'oscillatore di van der Pol\"\n");
    fprintf(gp, "set xlabel 'x1 (Tensione sul condensatore)'\n");
    fprintf(gp, "set ylabel 'x2 (Corrente sul condensatore)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key outside right top\n"); // Posizione della legenda (opzionale)

    std::string plot_cmd = "plot ";
    bool primo_plot = true;

    for (const auto &[chiave, evoluzione] : evoluzioni)
    {
        if (!primo_plot)
        {
            plot_cmd += ", ";
        }

        char title_buf[100];
        snprintf(title_buf, sizeof(title_buf), "x1_0 = %.2f, x2_0 = %.2f", chiave.first, chiave.second);

        plot_cmd += "'-' with lines title '" + std::string(title_buf) + "'";
        primo_plot = false;
    }
    fprintf(gp, "%s\n", plot_cmd.c_str());

    for (const auto &[chiave, evoluzione] : evoluzioni)
    {
        const auto &x1 = evoluzione.x1;
        const auto &x2 = evoluzione.x2;
        // Entrambi i vettori x1 e x2 hanno la stessa dimensione (N+1)
        for (size_t j = 0; j < x1.size(); ++j)
        {
            fprintf(gp, "%f %f\n", x1[j], x2[j]);
        }
        fprintf(gp, "e\n"); // Fine blocco dati per la curva corrente
    }

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    return 0;
}
