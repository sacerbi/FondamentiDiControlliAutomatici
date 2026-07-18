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
struct Params {
    double R;
    double C;
    double Vg;
    double w;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_circuitoRC(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    double u = p.Vg * sin(p.w * t); // ingresso sinusoidale
    double dx = (1/(p.R * p.C)) * (u - x[0]);
    return {dx};
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double R = 1, C = 1, Vg = 2, w = 3, x0 = 10;
    double dt = 0.001;
    double t0 = 0, tf = 20.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;
    std::vector<double> u(N + 1);
    for (int i = 0; i <= N; ++i) u[i] = Vg * sin(w * time[i]);

    // -----------------------------------------
    // Andamento della tensione sul condensatore nel circuito RC e della tensione sul resistore
    // -----------------------------------------
    double dummy = 0.0; // to match function signature
    auto y = rungeKutta4(f_circuitoRC, {x0}, t0, tf, dt, Params{R, C, Vg, w}, dummy);

    std::vector<double> Vc, Vr;
    for (auto& v : y) Vc.push_back(v[0]);
    for (int i = 0; i <= N; ++i) Vr.push_back(u[i] - Vc[i]);

    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE* gp = POPEN("gnuplot -persist", "w");
    if (!gp) {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return 1;
    }

    // ==========================================
    // FIGURA 2.4: Comportamento del circuito RC
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura 2.4'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Comportamento del circuito RC'\n");
    fprintf(gp, "set xlabel 'Tempo (t)'\n");
    fprintf(gp, "set ylabel 'Tensione (V)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    // Diciamo a Gnuplot di aspettarsi DUE serie di dati da stdin ('-' per Vc e '-' per Vr)
    fprintf(gp, "plot '-' with lines title 'V_c', '-' with lines title 'V_r'\n");

    // 1. Invio dei dati per la PRIMA curva (Vc)
    for (size_t i = 0; i < time.size(); ++i) {
        fprintf(gp, "%f %f\n", time[i], Vc[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva (Vc)

    // 2. Invio dei dati per la SECONDA curva (Vr)
    for (size_t i = 0; i < time.size(); ++i) {
        fprintf(gp, "%f %f\n", time[i], Vr[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la seconda curva (Vr)

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    std::vector<double> Vc_teorica(N + 1), Vr_teorica(N + 1);
    for (int i = 0; i <= N; ++i) {
        double gamma = atan(w * R * C);
        double mode = exp(-time[i]/(R*C));
        double modeGain = (Vg * w * R * C)/(1 + w*w*R*R*C*C);
        double amp = Vg / sqrt(1 + w*w*R*R*C*C);    
        Vc_teorica[i] = mode * x0 + modeGain * mode + amp * sin(w * time[i] - gamma);
        Vr_teorica[i] = -mode * x0 - modeGain * mode + amp * w * R * C * cos(w * time[i] - gamma);
    }

    // =================================================
    // FIGURA 2.4: Comportamento teorico del circuito RC
    // =================================================
    fprintf(gp, "set terminal qt 2 title 'Figura 2.4'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Comportamento del circuito RC'\n");
    fprintf(gp, "set xlabel 'Tempo (t)'\n");
    fprintf(gp, "set ylabel 'Tensione (V)'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    // Diciamo a Gnuplot di aspettarsi DUE serie di dati da stdin ('-' per Vc e '-' per Vr)
    fprintf(gp, "plot '-' with lines title 'V_c teorica', '-' with lines title 'V_r teorica'\n");

    // 1. Invio dei dati per la PRIMA curva (Vc)
    for (size_t i = 0; i < time.size(); ++i) {
        fprintf(gp, "%f %f\n", time[i], Vc_teorica[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva (Vc)

    // 2. Invio dei dati per la SECONDA curva (Vr)
    for (size_t i = 0; i < time.size(); ++i) {
        fprintf(gp, "%f %f\n", time[i], Vr_teorica[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la seconda curva (Vr)

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    double error_Vc = 0.0;
    double error_Vr = 0.0;
    for (int i = 0; i <= N; ++i) {
        error_Vc += std::pow(Vc[i] - Vc_teorica[i], 2);
        error_Vr += std::pow(Vr[i] - Vr_teorica[i], 2);
    }
    error_Vc = std::sqrt(error_Vc / (N + 1));
    error_Vr = std::sqrt(error_Vr / (N + 1));

    std::cout << "Errore RMS Vc: " << error_Vc << std::endl;
    std::cout << "Errore RMS Vr: " << error_Vr << std::endl;

    return 0;
}
