#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdio> // Per POPEN, PCLOSE, fprintf
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
// Parametri del sistema
// --------------------------
struct Params
{
    double u;
};
struct EmptyParams
{
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_sistema(double t, const std::vector<double> &x, const Params &p, const EmptyParams &p2)
{
    std::vector<double> dx(2);
    dx[0] = -x[0];
    dx[1] = x[1] * (x[0] - p.u);
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main()
{
    // Parametri di simulazione
    Params p{1};
    EmptyParams p2;
    double dt = 0.1;
    double ti = 0.0, tf = 25.0;

    // Vettori degli stati iniziali
    std::vector<double> x1_0 = {-3.0, -3.0, -3.0, -3.0, -3.0, -1.5, -1.5, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 3.0};
    std::vector<double> x2_0 = {-3.0, -1.0, 0.0, 1.0, 3.0, -3.0, 3.0, -3.0, 3.0, -3.0, 3.0, -3.0, 3.0, -3.0, 3.0, -3.0, -1.0, 0.0, 1.0, 3.0};

    // Struttura per memorizzare tutte le traiettorie in RAM prima di inviarle a Gnuplot
    std::vector<std::vector<std::vector<double>>> tutte_le_traiettorie;

    std::cout << "Calcolo delle traiettorie in corso...\n";
    // Simulazione e salvataggio in memoria
    for (size_t i = 0; i < x1_0.size(); i++)
    {
        std::vector<double> x0_vec = {x1_0[i], x2_0[i]};
        auto dx = rungeKutta4(f_sistema, x0_vec, ti, tf, dt, p, p2);
        tutte_le_traiettorie.push_back(dx);
    }
    std::cout << "Dati generati con successo. Avvio di Gnuplot...\n";

    // Avvio ed esecuzione di Gnuplot via Pipe
    FILE *gp = POPEN("gnuplot -persist", "w");
    if (!gp)
    {
        std::cerr << "Errore: impossibile avviare Gnuplot. Assicurati che sia installato e nel PATH di sistema." << std::endl;
        return 1;
    }

    // Impostazioni estetiche del grafico
    fprintf(gp, "set title 'Traiettorie e linee di livello (di Lyapunov)' font ',12'\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set size square\n");
    fprintf(gp, "set xrange [-3:3]\n");
    fprintf(gp, "set yrange [-3:3]\n");
    fprintf(gp, "set key off\n");

    std::cout << "Generazione curve di livello...\n";
    // Disegno le curve di livello
    std::vector<double> livelli_v = {0.125, 0.25, 2.0, 4.0};
    for (size_t i = 0; i < livelli_v.size(); ++i)
    {
        double raggio = std::sqrt(2.0 * livelli_v[i]);
        fprintf(gp, "set object %zu circle at 0,0 size %f fc rgb 'red' fs empty lw 1.5\n", i + 1, raggio);
    }

    std::vector<double> livelli_v2 = {-2.0, 2.0};
    for (size_t i = 0; i < livelli_v2.size(); ++i)
    {
        if (livelli_v2[i] > 0)
        { // <-- Aggiunto controllo per evitare sqrt di numeri negativi (NaN)
            double raggio = std::sqrt(2.0 * livelli_v2[i]);
            fprintf(gp, "set object %zu circle at 0,0 size %f fc rgb 'green' fs empty lw 1.5\n", i + 1 + livelli_v.size(), raggio);
        }
    }

    // ==========================================
    // COSTRUZIONE DEL COMANDO PLOT E INVIO DATI INLINE
    // ==========================================
    std::cout << "Generazione delle traiettorie...\n";
    fprintf(gp, "plot ");
    for (size_t i = 0; i < x1_0.size(); ++i)
    {
        // Prepariamo Gnuplot a ricevere x1_0.size() blocchi di dati inline ('-')
        fprintf(gp, "'-' with lines lc rgb 'blue' lw 1.2 notitle");
        if (i < x1_0.size() - 1)
        {
            fprintf(gp, ", ");
        }
    }
    fprintf(gp, "\n");

    // Invio sequenziale dei blocchi di dati calcolati
    for (size_t i = 0; i < tutte_le_traiettorie.size(); ++i)
    {
        for (const auto &x : tutte_le_traiettorie[i])
        {
            fprintf(gp, "%f %f\n", x[0], x[1]);
        }
        fprintf(gp, "e\n"); // 'e' indica a Gnuplot la fine della curva corrente
    }

    fflush(gp); // Forza l'invio dei comandi a Gnuplot
    PCLOSE(gp); // Chiude la pipe

    std::cout << "Grafico generato con successo!" << std::endl;
    return 0;
}