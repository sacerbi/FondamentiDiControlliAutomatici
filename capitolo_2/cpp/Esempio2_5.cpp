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
    double M;
    double Fm;
    double kt0;
    double alpha;
    double t0;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_sistemaCarrelloMolla(double t, const std::vector<double>& x, const Params& p, double& dummy) {
    double u = p.Fm;
    double Fe = p.kt0 * exp(-p.alpha*(t-p.t0)) * x[0];
    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = (u - Fe)/p.M;
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double Fm = 5, kt0 = 10, alpha = 2, s0 = 0, s0dot = 0, M = 10;
    double dt = 0.01;
    double t0 = 0, tf = 10.0;
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i) time[i] = t0 + i * dt;

    // -----------------------------------------
    // Andamento della tensione sul condensatore nel circuito RC e della tensione sul resistore
    // -----------------------------------------
    double dummy = 0.0; // to match function signature
    auto y = rungeKutta4(f_sistemaCarrelloMolla, {s0, s0dot}, t0, tf, dt, Params{M, Fm, kt0, alpha, t0}, dummy);

    std::vector<double> s, sdot, et;
    for (auto& v : y) {
        s.push_back(v[0]);
        sdot.push_back(v[1]);
    }
    for (int i = 0; i <= N; ++i){
        double Ue = (kt0 * exp(-alpha*(time[i]-t0)) * s[i]*s[i])/2;
        double Up = (M*sdot[i]*sdot[i])/2;
        et.push_back(Ue+Up);
    }

    // Apre una pipe verso Gnuplot (-persist mantiene la finestra aperta alla fine)
    FILE* gp = POPEN("gnuplot -persist", "w");
    if (!gp) {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return 1;
    }

    // ==========================================
    // FIGURA: Comportamento del sistema carrello-molla
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura'\n"); // Apre la finestra Qt 1
    fprintf(gp, "set title 'Comportamento del sistema carrello-molla'\n");
    fprintf(gp, "set xlabel 'Tempo (t)'\n");
    fprintf(gp, "set ylabel 'Spostamento / Velocità / Energia Totale'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set key top right\n"); // Posizione della legenda (opzionale)

    fprintf(gp, "plot '-' with lines title 's(t)', '-' with lines title 'sdot(t)', '-' with lines title 'E_T(t)'\n");

    // 1. Invio dei dati per la PRIMA curva
    for (size_t i = 0; i < time.size(); ++i) {
        fprintf(gp, "%f %f\n", time[i], s[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la prima curva

    // 2. Invio dei dati per la SECONDA curva
    for (size_t i = 0; i < time.size(); ++i) {
        fprintf(gp, "%f %f\n", time[i], sdot[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la seconda curva

    // 3. Invio dei dati per la TERZA curva
    for (size_t i = 0; i < time.size(); ++i) {
        fprintf(gp, "%f %f\n", time[i], et[i]);
    }
    fprintf(gp, "e\n"); // 'e' comunica la fine dei dati per la terza curva

    // Fondamentale: forza lo svuotamento del buffer per disegnare subito il grafico
    fflush(gp);

    return 0;
}
