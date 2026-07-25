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
    double dx = x[0] > 0 ? (p.ubar - 2 * (1 - exp(-x[0]))) : (p.ubar + x[0]);
    return {dx};
}

template <typename Func>
std::pair<std::vector<double>, std::vector<double>> evolviSistema(Func sistema, double ti, double tf, double dt, double x0, const Params &p)
{
    int M = static_cast<int>((tf - ti) / dt);
    std::vector<double> t(M + 1);
    for (int i = 0; i <= M; ++i)
        t[i] = ti + i * dt;
    double dummy = 0.0; // to match function signature
    auto y = rungeKutta4(sistema, {x0}, ti, tf, dt, p, dummy);
    std::vector<double> y_vals;
    for (auto &v : y)
        y_vals.push_back(v[0]);
    return {t, y_vals};
}

void esercizio2_7_1(std::vector<double> x)
{
    // Dati
    double u1 = -1.0, u2 = 0.0, u3 = 1.0;
    std::vector<double> xp1(x.size());
    std::vector<double> xp2(x.size());
    std::vector<double> xp3(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
        xp1[i] = u1 * (x[i] * x[i] * x[i]);
        xp2[i] = u2 * (x[i] * x[i] * x[i]);
        xp3[i] = u3 * (x[i] * x[i] * x[i]);
    }
    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return;
    }

    // ==========================================
    // FIGURA 2.7.1
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura 2.7.1'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Esercizio 2.7.1'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'dx/dt'\n");
    fprintf(gp, "set xrange [%f:%f]\n", x.front(), x.back());
    fprintf(gp, "set yrange [%f:%f]\n", -1.0, 2.0);
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'ubar = -1.0', '-' with lines title 'ubar = 0.0', '-' with lines title 'ubar = 1.0'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xp1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xp2[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
                        // 3. Invio dei dati per la TERZA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xp3[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
}

void esercizio2_7_2(std::vector<double> x)
{
    // Dati
    double u1 = -1.0, u2 = 0.0, u3 = 1.0, u4 = 2.1;
    std::vector<double> xp1(x.size());
    std::vector<double> xp2(x.size());
    std::vector<double> xp3(x.size());
    std::vector<double> xp4(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
        if (x[i] <= 0)
        {
            xp1[i] = u1 + x[i];
            xp2[i] = u2 + x[i];
            xp3[i] = u3 + x[i];
            xp4[i] = u4 + x[i];
        }
        else
        {
            xp1[i] = u1 - 2 * (1 - exp(-x[i]));
            xp2[i] = u2 - 2 * (1 - exp(-x[i]));
            xp3[i] = u3 - 2 * (1 - exp(-x[i]));
            xp4[i] = u4 - 2 * (1 - exp(-x[i]));
        }
    }
    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return;
    }

    // ==========================================
    // FIGURA 2.7.2
    // ==========================================
    fprintf(gp, "set terminal qt 2 title 'Figura 2.7.2'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Esercizio 2.7.2 - Equilibrio'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'dx/dt'\n");
    fprintf(gp, "set xrange [%f:%f]\n", x.front(), x.back());
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'ubar = -1.0', '-' with lines title 'ubar = 0.0', '-' with lines title 'ubar = 1.0', '-' with lines title 'ubar = 2.1', 0 with lines title 'x(t) = 0'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xp1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xp2[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // 3. Invio dei dati per la TERZA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xp3[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // 4. Invio dei dati per la QUARTA curva
    for (size_t i = 0; i < x.size(); ++i)
    {
        fprintf(gp, "%f %f\n", x[i], xp4[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));

    auto [t1, xdot1] = evolviSistema(f_sistema, 0.0, 10.0, 0.001, -0.9, Params{1.0});
    auto [t2, xdot2] = evolviSistema(f_sistema, 0.0, 10.0, 0.001, 3.0, Params{1.0});
    auto xbar = -log(0.5);

    // ==========================================
    // FIGURA 2.7.2
    // ==========================================
    fprintf(gp, "set terminal qt 3 title 'Figura 2.7.2'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Esercizio 2.7.2 - Stabilita'\n");
    fprintf(gp, "set xlabel 't'\n");
    fprintf(gp, "set ylabel 'x(t)'\n");
    fprintf(gp, "set xrange [%f:%f]\n", t1.front(), t1.back());
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 'u(t) = 1, x(0) = -0.9', '-' with lines title 'u(t) = 1, x(0) = 3.0', %f with lines title 'x = %f'\n", xbar, xbar);

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < t1.size(); ++i)
    {
        fprintf(gp, "%f %f\n", t1[i], xdot1[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < t1.size(); ++i)
    {
        fprintf(gp, "%f %f\n", t1[i], xdot2[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva
    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
}

void esercizio2_7_3(std::vector<double> x)
{
    std::map<std::pair<double, double>, std::vector<double>> risultati;
    // Caso ubar = 0, diversi valori di alpha
    double ubar = 0;
    std::vector<double> alphas = {-1.0, 0.0, 1.0};
    for (auto alpha : alphas)
    {
        std::vector<double> xdot(x.size());
        for (int i = 0; i < xdot.size(); i++)
        {
            xdot[i] = alpha * x[i] - ubar * ubar * x[i] * x[i] * x[i];
        }
        risultati[{ubar, alpha}] = xdot;
    }
    // Caso ubar != 0, alpha diversi
    ubar = 1.0;
    for (auto alpha : alphas)
    {
        std::vector<double> xdot(x.size());
        for (int i = 0; i < xdot.size(); i++)
        {
            xdot[i] = alpha * x[i] - ubar * ubar * x[i] * x[i] * x[i];
        }
        risultati[{ubar, alpha}] = xdot;
    }

    // Plot dei risultati
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return;
    }
    fprintf(gp, "set terminal qt 4 title 'Figura 2.7.3'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Esercizio 2.7.3'\n");
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'dx/dt'\n");
    fprintf(gp, "set xrange [%f:%f]\n", x.front(), x.back());
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n");

    std::string plot_cmd = "plot ";
    bool primo_plot = true;

    for (const auto &[chiave, evoluzione] : risultati)
    {
        if (!primo_plot)
        {
            plot_cmd += ", ";
        }
        // Aggiungiamo '-' per indicare che i dati arriveranno da standard input
        plot_cmd += "'-' with lines title 'u = " + std::to_string(chiave.first) +
                    ", alpha = " + std::to_string(chiave.second) + "'";
        primo_plot = false;
    }

    // Aggiungiamo la retta y=0 alla fine del comando
    // Gnuplot può plottare la costante '0' direttamente senza bisogno di passargli i dati (x, 0)
    if (!primo_plot)
    {
        plot_cmd += ", ";
    }
    plot_cmd += "0 with lines title 'y = 0' linecolor 'black' linetype 1";

    // 2. Inviamo il comando plot completo a Gnuplot
    fprintf(gp, "%s\n", plot_cmd.c_str());

    // 3. Inviamo sequenzialmente i dati per ogni curva definita con '-'
    for (const auto &[chiave, evoluzione] : risultati)
    {
        for (size_t i = 0; i < x.size() && i < evoluzione.size(); ++i)
        {
            fprintf(gp, "%f %f\n", x[i], evoluzione[i]);
        }
        fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la curva corrente
    }
    fflush(gp);
}

// --------------------------
// MAIN
// --------------------------
int main()
{
    // Dati
    double u1 = 0.125, u2 = 0.25;
    double dx = 0.01;
    double x0 = -3, xf = 3;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i)
        x[i] = x0 + i * dx;

    // Punto 1
    esercizio2_7_1(x);

    // Punto 2
    esercizio2_7_2(x);

    // Punto 3
    esercizio2_7_3(x);

    return 0;
}
