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
    double g, k, L, M;
    std::vector<double> u;
    double dt;
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_pendolo(double t, const std::vector<double> &x, const Params &p, double &dummy)
{
    int idx = static_cast<int>(t / p.dt);
    double ut = idx < p.u.size() ? p.u[idx] : p.u.back();

    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = -(p.g / p.L) * sin(x[0]) - (p.k / (p.M * p.L * p.L)) * x[1] + 1 / (p.M * p.L * p.L) * ut;
    return dx;
}

int main()
{
    // Dati
    double M = 10.0, k = 2.0, g = 9.81, L = 0.5;
    double dt = 0.001;
    double ti = 0.0, tf = 10.0;
    std::vector<double> x0 = {M_PI / 2, 0.0};
    int N = static_cast<int>((tf - ti) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i)
        time[i] = ti + i * dt;
    std::vector<double> u(N + 1);
    for (int i = 0; i <= N; ++i)
        u[i] = 0.0;
    Params p1 = {g, k, L, M, u, dt};

    // Simulazione del sistema
    double dummy = 0.0;
    auto dx1 = rungeKutta4(f_pendolo, x0, ti, tf, dt, p1, dummy);

    std::vector<double> x1, x2, y;
    for (auto &x : dx1)
    {
        x1.push_back(x[0]);
        x2.push_back(x[1]);
        double x2_val = (1 / 2) * M * L * L * x[1] * x[1];
        double x1_val = -M * g * L * cos(x[0]);
        y.push_back(x2_val + x1_val);
    }

    // Plot Caso 1
    // --- 1. SALVIAMO I DATI IN UN FILE ---
    std::ofstream dataFile("esempio3_18.dat");
    for (size_t i = 0; i < time.size(); i++)
    {
        // Scriviamo in colonne: tempo | E | pos | vel
        dataFile << time[i] << " " << y[i] << " " << x1[i] << " " << x2[i] << "\n";
    }
    dataFile.close();

    // --- 2. CREIAMO LE ISTRUZIONI PER IL GRAFICO ---
    std::ofstream scriptFile("plot.gp");
    scriptFile << "set title 'Comportamento del sistema pendolo'\n";
    scriptFile << "set xlabel 'Tempo (t)'\n";
    scriptFile << "set ylabel 'Posizione/Velocità/Energia'\n";
    scriptFile << "set grid\n";
    // Plottiamo entrambe le curve nello stesso grafico per confrontarle
    scriptFile << "plot 'esempio3_18.dat' using 1:2 with lines lw 2 title \"E\", \\\n";
    scriptFile << "     'esempio3_18.dat' using 1:3 with lines lw 2 title \"pos\", \\\n"; // <- virgola!
    scriptFile << "     'esempio3_18.dat' using 1:4 with lines lw 2 title \"vel\"\n";     // <- colonna 4, non 3
    scriptFile.close();

    // --- 3. LANCIAMO GNUPLOT ---
    std::cout << "Apertura del grafico in corso..." << std::endl;
    // Il flag -p (persist) dice a gnuplot di lasciare aperta la finestra del grafico
    system("gnuplot -p plot.gp");

    std::cout << "Premi INVIO per continuare..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::remove("esempio3_18.dat");
    std::remove("plot.gp");

    return 0;
}