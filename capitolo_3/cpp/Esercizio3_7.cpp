#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <cmath>
#include "myDiffEquation.h"
#include <fstream> // <--- AGGIUNTO per salvare i file
#include <cstdlib> // <--- AGGIUNTO per lanciare Gnuplot

// --------------------------
// Parametri del sistema
// --------------------------
struct Params
{
    double A, w, ubar;
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_sistema(double t, const std::vector<double> &x, const Params &p, double &dummy)
{
    double ut = p.ubar + p.A * cos(p.w * t);
    std::vector<double> dx(1);
    dx[0] = x[0] * x[0] - ut * x[0] - 2 * ut;
    return dx;
}

int main()
{
    // Dati
    double A = 0.1, w = 2.0, ubar = 1.0;
    double dt = 0.001;
    double ti = 0.0, tf = 10.0;
    std::vector<double> x0 = {-1.0};
    int N = static_cast<int>((tf - ti) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i)
        time[i] = ti + i * dt;
    Params p1 = {A, w, ubar};

    // Simulazione del sistema
    double dummy = 0.0;
    auto dx = rungeKutta4(f_sistema, x0, ti, tf, dt, p1, dummy);
    std::vector<double> y;
    for (std::size_t i = 0; i < dx.size(); i++)
    {
        auto &x = dx[i];
        double yt = pow(x[0], 3.0) + pow(ubar + A * cos(w * time[i]), 3.0);
        y.push_back(yt);
    }

    // Sistema linearizzato
    std::vector<double> ulin(N + 1);
    for (int i = 0; i <= N; i++)
        ulin[i] = A * std::cos(w * time[i]);
    std::vector<std::vector<double>> mA = {{-3.0}};
    std::vector<std::vector<double>> mB = {{-1.0}};
    std::vector<std::vector<double>> mC = {{3.0}};
    std::vector<std::vector<double>> mD = {{3.0}};

    std::vector<double> yl = lsim(mA, mB, mC, mD, ulin, time);

    // Plot Caso 1
    // --- 1. SALVIAMO I DATI IN UN FILE ---
    std::ofstream dataFile("esercizio3_7.dat");
    for (size_t i = 0; i < time.size(); i++)
    {
        // Scriviamo in colonne: tempo | y | y linearizzato
        dataFile << time[i] << " " << y[i] << " " << yl[i] << "\n";
    }
    dataFile.close();

    // --- 2. CREIAMO LE ISTRUZIONI PER IL GRAFICO ---
    std::ofstream scriptFile("plot.gp");
    scriptFile << "set title 'Comportamento del sistema'\n";
    scriptFile << "set xlabel 'Tempo (t)'\n";
    scriptFile << "set ylabel 'Risposta'\n";
    scriptFile << "set grid\n";
    // Plottiamo entrambe le curve nello stesso grafico per confrontarle
    scriptFile << "plot 'esercizio3_7.dat' using 1:2 with lines lw 2 title \"y\", \\\n";
    scriptFile << "     'esercizio3_7.dat' using 1:3 with lines lw 2 title \"y lin\"\n"; // <- colonna 4, non 3
    scriptFile.close();

    // --- 3. LANCIAMO GNUPLOT ---
    std::cout << "Apertura del grafico in corso..." << std::endl;
    // Il flag -p (persist) dice a gnuplot di lasciare aperta la finestra del grafico
    system("gnuplot -p plot.gp");

    std::cout << "Premi INVIO per continuare..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::remove("esercizio3_7.dat");
    std::remove("plot.gp");

    return 0;
}