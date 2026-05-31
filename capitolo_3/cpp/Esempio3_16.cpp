#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream> // <--- AGGIUNTO per salvare i file
#include <cstdlib> // <--- AGGIUNTO per lanciare Gnuplot

// ----------------------------------
// Funzione ausiliaria per verificare che una riga sia nulla
// ----------------------------------
bool isRowZero(const std::vector<double> &v, double epsilon = 0.0001)
{
    for (auto val : v)
    {
        if (std::abs(val) > epsilon)
            return false;
    }
    return true;
}

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<std::vector<double>> routhTable(const std::vector<double> &poly, double epsilon = 0.0001)
{
    int nCoeff = poly.size();
    std::cout << "Routh Table method for " << nCoeff - 1 << " degree polynomial. Using epsilon = " << epsilon << "\n";

    // Creazione della tabella vuota
    std::vector<std::vector<double>> rt(nCoeff, std::vector<double>(static_cast<int>(ceil(nCoeff / 2.0)), 0.0));

    // Popolo le prime due righe
    for (int i = 0; i < nCoeff; i++)
    {
        int r = i % 2;                         // Decido se prima o seconda riga a seconda che l'esponente sia pari o dispari
        int c = static_cast<int>(ceil(i / 2)); // Procedo con le colonne
        rt[r][c] = poly[i];
    }

    int rowsToDo = nCoeff - 2;
    std::vector<int> maxIndexColumn(rowsToDo, 0);
    for (int i = 0; i < rowsToDo; i++)
        maxIndexColumn[rowsToDo - i - 1] = static_cast<int>(ceil(i / 2)) + 1;

    for (int i = 2; i < nCoeff; i++)
    {
        if (isRowZero(rt[i - 1]))
        {
            std::cout << "\nCaso speciale riga nulla\n";
            int orderAuxPoly = nCoeff - i + 1;                          // Ordine del polinomio ausiliario
            int nCoeffAux = static_cast<int>(ceil(orderAuxPoly / 2.0)); // Num. coeff. polinomio ausiliario
            std::vector<double> auxPoly(nCoeffAux);                     // Prendo il polinomio ausiliario da usare
            for (int p = 0; p < nCoeffAux; p++)
            {
                auxPoly[p] = rt[i - 2][p];
            }
            std::vector<int> auxPow; // Il polinomio ha potenze pari
            for (int k = orderAuxPoly; k >= 0; k -= 2)
            {
                auxPow.push_back(k);
            }
            // Assegno la derivata
            for (int k = 0; k < nCoeffAux; k++)
            {
                rt[i - 1][k] = auxPoly[k] * auxPow[k];
            }
        }
        else if (std::abs(rt[i - 1][0]) < 0.0001)
        {
            std::cout << "\nCaso speciale primo termine nullo\n";
            bool allPositive = true;
            for (int p = 0; p < i - 1; p++)
            {
                if (!allPositive)
                    break;
                allPositive &= rt[p][0] >= 0;
            }
            rt[i - 1][0] = allPositive ? epsilon : -epsilon;
        }
        for (int j = 0; j < maxIndexColumn[i - 2]; j++)
        {
            double h1 = rt[i - 2][0];
            double k1 = rt[i - 1][0];
            double hi1 = rt[i - 2][j + 1];
            double ki1 = rt[i - 1][j + 1];
            rt[i][j] = hi1 - ((h1 * ki1) / k1);
        }
    }

    return rt;
}

bool isStabile(const std::vector<std::vector<double>> &rt)
{
    bool tuttiPositivi = true;
    bool tuttiNegativi = true;
    for (const auto &riga : rt)
    {
        if (riga[0] <= 0)
            tuttiPositivi = false;
        if (riga[0] >= 0)
            tuttiNegativi = false;
    }
    return tuttiPositivi || tuttiNegativi;
}

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
            std::vector<double> poli = {1.0, 2.0 + b, 1.0 + 2.0 * b, a + b};

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