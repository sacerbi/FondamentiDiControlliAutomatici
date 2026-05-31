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
    double k, h, M;
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_massSpringDamper(double t, const std::vector<double> &x, const Params &p, double &dummy)
{
    std::vector<double> dx(2);
    dx[0] = x[1];
    dx[1] = -(p.k / p.M) * x[0] - (p.h / p.M) * x[1] + (1 / p.M) * 0.0; //  Nessuna forza esterna
    return dx;
}

std::vector<double> f_caso1(const std::vector<double> &t, const std::vector<double> &x0, const Params &p)
{
    double disc = (p.h * p.h) / (4 * p.M * p.M) - p.k / p.M;
    double s1 = -p.h / (2 * p.M) + std::sqrt(disc);
    double s2 = -p.h / (2 * p.M) - std::sqrt(disc);
    std::vector<double> y_teorica;
    for (int i = 0; i < t.size(); i++)
    {
        double val = (1.0 / (s2 - s1)) * ((s2 * x0[0] - x0[1]) * exp(s1 * t[i]) - (s1 * x0[0] - x0[1]) * exp(s2 * t[i]));
        y_teorica.push_back(val);
    }
    return y_teorica;
}

std::vector<double> f_caso2(const std::vector<double> &t, const std::vector<double> &x0, const Params &p)
{
    double sigma = -p.h / (2 * p.M);
    double omega = std::sqrt(p.k / p.M - (p.h * p.h) / (4 * p.M * p.M));
    double delta = std::sqrt(sigma * sigma + omega * omega);
    double gamma = std::acos(sigma / delta);
    std::vector<double> y_teorica;
    for (int i = 0; i < t.size(); i++)
    {
        double val = exp(sigma * t[i]) * (-(delta / omega) * std::sin(omega * t[i] - gamma) * x0[0] + (1 / omega) * std::sin(omega * t[i]) * x0[1]);
        y_teorica.push_back(val);
    }
    return y_teorica;
}

std::vector<double> f_caso3(const std::vector<double> &t, const std::vector<double> &x0, const Params &p)
{
    double s0 = -p.h / (2 * p.M);
    std::vector<double> y_teorica;
    for (int i = 0; i < t.size(); i++)
    {
        double val = exp(s0 * t[i]) * (x0[0] - (s0 * x0[0] - x0[1]) * t[i]);
        y_teorica.push_back(val);
    }
    return y_teorica;
}

double rms(const std::vector<double> &s1, const std::vector<double> &s2, double n)
{
    for (int i = 0; i < s1.size(); ++i)
    {
        n += std::pow(s1[i] - s2[i], 2);
    }
    n = std::sqrt(n / (s1.size() + 1));
    return n;
}

int main()
{
    // Configurazione Caso 1: h^2 > 4Mk
    Params p1 = {1.0, 3.0, 2.0};
    double dt = 0.001;
    double ti = 0.0, tf = 50.0;
    std::vector<double> x0 = {2.0, 1.0};
    int N = static_cast<int>((tf - ti) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i)
        time[i] = ti + i * dt;

    // Simulazione del sistema
    double dummy = 0.0;
    auto dx1 = rungeKutta4(f_massSpringDamper, x0, ti, tf, dt, p1, dummy);

    std::vector<double> y1;
    for (auto &x : dx1)
        y1.push_back(x[0]);

    // Confronto del sistema con la soluzione teorica
    std::vector<double> y1_teorica = f_caso1(time, x0, p1);

    // Plot Caso 1
    // --- 1. SALVIAMO I DATI IN UN FILE ---
    std::ofstream dataFile("esempio3_5.dat");
    for (size_t i = 0; i < time.size(); i++)
    {
        // Scriviamo in colonne: tempo | y1 | y1 teorica
        dataFile << time[i] << " " << y1[i] << " " << y1_teorica[i] << "\n";
    }
    dataFile.close();

    // --- 2. CREIAMO LE ISTRUZIONI PER IL GRAFICO ---
    std::ofstream scriptFile1("plot.gp");
    scriptFile1 << "set title 'Comportamento del sistema k,h,M'\n";
    scriptFile1 << "set xlabel 'Tempo (t)'\n";
    scriptFile1 << "set ylabel 'Posizione (m)'\n";
    scriptFile1 << "set grid\n";
    // Plottiamo entrambe le curve nello stesso grafico per confrontarle
    scriptFile1 << "plot 'esempio3_5.dat' using 1:2 with lines lw 2 title \"y1\", \\\n";
    scriptFile1 << "     'esempio3_5.dat' using 1:3 with lines lw 2 title \"y1 teorica\"\n";
    scriptFile1.close();

    // --- 3. LANCIAMO GNUPLOT ---
    std::cout << "Apertura del grafico in corso..." << std::endl;
    system("gnuplot -p plot.gp");

    std::cout << "Premi INVIO per continuare..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::remove("esempio3_6.dat");
    std::remove("plot.gp");

    // Calcolo dell'errore
    double error1 = rms(y1, y1_teorica, 0.0);
    std::cout << "Errore RMS y: " << error1 << std::endl;

    // Configurazione Caso 2: h^2 < 4Mk
    Params p2 = {1.0, 1.0, 2.0};
    auto dx2 = rungeKutta4(f_massSpringDamper, x0, ti, tf, dt, p2, dummy);

    std::vector<double> y2;
    for (auto &x : dx2)
        y2.push_back(x[0]);

    // Confronto del sistema con la soluzione teorica
    std::vector<double> y2_teorica = f_caso2(time, x0, p2);

    // Plot Caso 2
    // --- 1. SALVIAMO I DATI IN UN FILE ---
    std::ofstream dataFile2("esempio3_5.dat");
    for (size_t i = 0; i < time.size(); i++)
    {
        // Scriviamo in colonne: tempo | y2 | y2 teorica
        dataFile2 << time[i] << " " << y2[i] << " " << y2_teorica[i] << "\n";
    }
    dataFile2.close();

    // --- 2. CREIAMO LE ISTRUZIONI PER IL GRAFICO ---
    std::ofstream scriptFile2("plot.gp");
    scriptFile2 << "set title 'Comportamento del sistema k,h,M'\n";
    scriptFile2 << "set xlabel 'Tempo (t)'\n";
    scriptFile2 << "set ylabel 'Posizione (m)'\n";
    scriptFile2 << "set grid\n";
    // Plottiamo entrambe le curve nello stesso grafico per confrontarle
    scriptFile2 << "plot 'esempio3_5.dat' using 1:2 with lines lw 2 title \"y2\", \\\n";
    scriptFile2 << "     'esempio3_5.dat' using 1:3 with lines lw 2 title \"y2 teorica\"\n";
    scriptFile2.close();

    // --- 3. LANCIAMO GNUPLOT ---
    std::cout << "Apertura del grafico in corso..." << std::endl;
    system("gnuplot -p plot.gp");

    std::cout << "Premi INVIO per continuare..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::remove("esempio3_6.dat");
    std::remove("plot.gp");

    // Calcolo dell'errore
    double error2 = rms(y2, y2_teorica, 0.0);
    std::cout << "Errore RMS y: " << error2 << std::endl;

    // Configurazione Caso 1: h^2 == 4Mk
    Params p3 = {2.0, 4.0, 2.0};
    auto dx3 = rungeKutta4(f_massSpringDamper, x0, ti, tf, dt, p3, dummy);

    std::vector<double> y3;
    for (auto &x : dx3)
        y3.push_back(x[0]);

    // Confronto del sistema con la soluzione teorica
    std::vector<double> y3_teorica = f_caso3(time, x0, p3);

    // Plot Caso 3
    // --- 1. SALVIAMO I DATI IN UN FILE ---
    std::ofstream dataFile3("esempio3_5.dat");
    for (size_t i = 0; i < time.size(); i++)
    {
        // Scriviamo in colonne: tempo | y1 | y1 teorica
        dataFile3 << time[i] << " " << y3[i] << " " << y3_teorica[i] << "\n";
    }
    dataFile3.close();

    // --- 2. CREIAMO LE ISTRUZIONI PER IL GRAFICO ---
    std::ofstream scriptFile3("plot.gp");
    scriptFile3 << "set title 'Comportamento del sistema k,h,M'\n";
    scriptFile3 << "set xlabel 'Tempo (t)'\n";
    scriptFile3 << "set ylabel 'Posizione (m)'\n";
    scriptFile3 << "set grid\n";
    // Plottiamo entrambe le curve nello stesso grafico per confrontarle
    scriptFile3 << "plot 'esempio3_5.dat' using 1:2 with lines lw 2 title \"y3\", \\\n";
    scriptFile3 << "     'esempio3_5.dat' using 1:3 with lines lw 2 title \"y3 teorica\"\n";
    scriptFile3.close();

    // --- 3. LANCIAMO GNUPLOT ---
    std::cout << "Apertura del grafico in corso..." << std::endl;
    system("gnuplot -p plot.gp");

    std::cout << "Premi INVIO per continuare..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::remove("esempio3_6.dat");
    std::remove("plot.gp");
    // Calcolo dell'errore
    double error3 = rms(y3, y3_teorica, 0.0);
    std::cout << "Errore RMS y: " << error3 << std::endl;

    return 0;
}