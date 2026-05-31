#include <iostream>
#include <iomanip>
#include "myDiffEquation.h"
#include <fstream> // <--- AGGIUNTO per salvare i file
#include <cstdlib> // <--- AGGIUNTO per lanciare Gnuplot

void printTable(const std::vector<std::vector<double>> &t)
{
    std::cout << "Table t:" << std::endl;
    for (auto v : t)
    {
        for (auto x : v)
        {
            std::cout << x << "\t";
        }
        std::cout << std::endl;
    }
}

int main()
{
    std::vector<double> x_stabile, y_stabile;
    std::vector<double> x_instabile, y_instabile;

    double start = -2.0, end = 4.0, step = 0.1;

    for (double a = start; a <= end; a += step)
    {
        for (double b = start; b <= end; b += step)
        {
            std::vector<double> poli = {1.0, a, b, 1.0};

            auto tabella = routhTable(poli);
            if (isStabile(tabella))
            {
                x_stabile.push_back(a);
                y_stabile.push_back(b);
            }
            else
            {
                x_instabile.push_back(a);
                y_instabile.push_back(b);
            }
        }
    }

    // --- 1. SALVIAMO I DATI IN UN FILE ---
    std::ofstream fStabile("stabile.dat");
    for (size_t i = 0; i < x_stabile.size(); i++)
        fStabile << x_stabile[i] << " " << y_stabile[i] << "\n";
    fStabile.close();

    std::ofstream fInstabile("instabile.dat");
    for (size_t i = 0; i < x_instabile.size(); i++)
        fInstabile << x_instabile[i] << " " << y_instabile[i] << "\n";
    fInstabile.close();

    // --- 2. CREIAMO LE ISTRUZIONI PER IL GRAFICO ---
    std::ofstream scriptFile("plot.gp");
    scriptFile << "set title 'Regione di Stabilita'\n";
    scriptFile << "set xlabel 'alpha'\n";
    scriptFile << "set ylabel 'beta'\n";
    scriptFile << "set key top right box opaque\n";
    // Plottiamo entrambe le curve nello stesso grafico per confrontarle
    scriptFile << "plot 'stabile.dat'   with points pt 1 ps 0.5 lc rgb 'red'  title 'Stabile',\\\n";
    scriptFile << "     'instabile.dat' with points pt 2 ps 0.5 lc rgb 'blue' title 'Instabile'\n";
    scriptFile.close();

    // --- 3. LANCIAMO GNUPLOT ---
    std::cout << "Apertura del grafico in corso..." << std::endl;
    // Il flag -p (persist) dice a gnuplot di lasciare aperta la finestra del grafico
    system("gnuplot -p plot.gp");

    std::cout << "Premi INVIO per continuare..." << std::endl;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::remove("stabile.dat");
    std::remove("instabile.dat");
    std::remove("plot.gp");

    return 0;
}