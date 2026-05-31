#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include "myDiffEquation.h"
#include <fstream> // <--- AGGIUNTO per salvare i file
#include <cstdlib> // <--- AGGIUNTO per lanciare Gnuplot

// --------------------------
// Parametri del sistema
// --------------------------
struct Params
{
    std::array<std::array<double, 2>, 2> A{};
    std::array<std::array<double, 1>, 2> B{};
    std::array<double, 2> C{};
    double D;
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_system(double t, const std::vector<double> &x, const Params &p, const std::function<double(double)> &u_func)
{
    std::vector<double> dx(2);
    double u = u_func(t);
    dx[0] = p.A[0][0] * x[0] + p.A[0][1] * x[1] + p.B[0][0] * u;
    dx[1] = p.A[1][0] * x[0] + p.A[1][1] * x[1] + p.B[1][0] * u;
    return dx;
}

bool isStable(const std::array<std::array<double, 2>, 2> &A)
{
    auto poles = eigenvalues2x2(A);
    std::cout << "Poli del sistema:\n";
    for (int i = 0; i < 2; ++i)
    {
        std::cout << "  p" << (i + 1) << " = "
                  << poles[i].real();
        if (poles[i].imag() != 0.0)
            std::cout << " + " << poles[i].imag() << "j";
        std::cout << "\n";
    }
    // Stabile se tutte le parti reali < 0
    return poles[0].real() < 0.0 && poles[1].real() < 0.0;
}

double staticGain(const Params &p)
{
    const auto &A = p.A;
    const auto &B = p.B;
    const auto &C = p.C;

    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];

    // Inversa di A (matrice 2x2)
    // inv(A) = (1/det) * [[ A[1][1], -A[0][1] ],
    //                     [-A[1][0],  A[0][0] ]]

    // Calcola (-A)^-1 * B  →  equivalente a  -(A^-1 * B)
    // invA_B = A^-1 * B
    double invA_B0 = (A[1][1] * B[0][0] - A[0][1] * B[1][0]) / det;
    double invA_B1 = (-A[1][0] * B[0][0] + A[0][0] * B[1][0]) / det;

    // K = C * (-A)^-1 * B + D  =  -C * A^-1 * B + D
    double K = -(C[0] * invA_B0 + C[1] * invA_B1) + p.D;
    return K;
}

int main()
{
    Params p = {
        {{{-1, 0}, {0, -2}}},
        {{{5}, {4}}},
        {1, 3},
        0};
    double dt = 0.001;
    double ti = 0.0, tf = 10.0;
    std::vector<double> x0 = {0.0, 0.0};
    int N = static_cast<int>((tf - ti) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        time[i] = ti + i * dt;
    }

    // Impulso: valore alto solo per t < dt, zero altrove
    auto u_imp = [dt](double t)
    { return t < dt ? 1.0 / dt : 0.0; };

    // Rampa: u(t) = t
    auto u_ramp = [](double t)
    { return t; };

    // Verifica di stabilità
    bool stabile = isStable(p.A);
    if (stabile)
        std::cout << "Sistema stabile" << std::endl;
    else
        std::cout << "Sistema NON stabile" << std::endl;

    // Simulazione della risposta all'impulso
    auto dx1 = rungeKutta4(f_system, x0, ti, tf, dt, p, u_imp);

    std::vector<double> y_imp;
    y_imp.reserve(dx1.size());
    for (std::size_t i = 0; i < dx1.size(); ++i)
    {
        double u = u_imp(time[i]);
        double y = p.C[0] * dx1[i][0] + p.C[1] * dx1[i][1] + p.D * u;
        y_imp.push_back(y);
    }

    // Simulazione della risposta alla rampa
    auto dx2 = rungeKutta4(f_system, x0, ti, tf, dt, p, u_ramp);

    std::vector<double> y_ramp;
    y_ramp.reserve(dx2.size());
    for (std::size_t i = 0; i < dx2.size(); ++i)
    {
        double u = u_ramp(time[i]);
        double y = p.C[0] * dx2[i][0] + p.C[1] * dx2[i][1] + p.D * u;
        y_ramp.push_back(y);
    }

    // --- 1. SALVIAMO I DATI IN UN FILE ---
    std::ofstream dataFile("esercizio3_4.dat");
    for (size_t i = 0; i < time.size(); i++)
    {
        // Scriviamo in colonne: tempo | impulso | rampa
        dataFile << time[i] << " " << y_imp[i] << " " << y_ramp[i] << "\n";
    }
    dataFile.close();

    // --- 2. CREIAMO LE ISTRUZIONI PER IL GRAFICO ---
    std::ofstream scriptFile("plot.gp");
    scriptFile << "set title 'Analisi del Sistema'\n";
    scriptFile << "set xlabel 'Tempo (t)'\n";
    scriptFile << "set ylabel 'Ampiezza'\n";
    scriptFile << "set grid\n";
    // Plottiamo entrambe le curve nello stesso grafico per confrontarle
    scriptFile << "plot 'esercizio3_4.dat' using 1:2 with lines lw 2 title \"Impulso\", \\\n";
    scriptFile << "     'esercizio3_4.dat' using 1:3 with lines lw 2 title \"Rampa\"\n";
    scriptFile.close();

    // --- 3. LANCIAMO GNUPLOT ---
    std::cout << "Apertura del grafico in corso..." << std::endl;
    // Il flag -p (persist) dice a gnuplot di lasciare aperta la finestra del grafico
    system("gnuplot -p plot.gp");

    // calcolo del guadagno statico
    auto K = staticGain(p);
    std::cout << "Guadagno statico: " << K << std::endl;

    std::cout << "Premi INVIO per continuare..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::remove("esercizio3_4.dat");
    std::remove("plot.gp");

    return 0;
}