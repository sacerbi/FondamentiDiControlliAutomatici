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
    dx[0] = -x[0] * (x[0] + p.u) + x[1] * x[1];
    dx[1] = -(x[0] + 2 * p.u) * x[1];
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
    double ti = 0.0, tf = 20.0, tf_diverge = 4.0;

    // Vettori degli stati iniziali
    std::vector<double> x1_0 = {-0.99, -2.0, -2.0, -1.5, -2.0, -2.0, -1.5, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0};
    std::vector<double> x2_0 = {0.0, 1.48, 1.8, 2.0, -1.48, -1.8, -2.0, 2.0, -2.0, 2.0, -2.0, 2.0, -2.0, 1.5, -1.5};
    std::vector<double> x1_0_diverge = {-2.0, -2.0, -2.0, -2.0, -1.01};
    std::vector<double> x2_0_diverge = {-1.45, 1.45, -1.42, 1.42, 0.0};

    std::vector<std::vector<std::vector<double>>> tutte_le_traiettorie;

    std::cout << "Calcolo delle traiettorie in corso...\n";
    // Simulazione e salvataggio in memoria
    for (size_t i = 0; i < x1_0.size(); i++)
    {
        std::vector<double> x0_vec = {x1_0[i], x2_0[i]};
        auto dx = rungeKutta4(f_sistema, x0_vec, ti, tf, dt, p, p2);
        tutte_le_traiettorie.push_back(dx);
    }
    for (size_t i = 0; i < x1_0_diverge.size(); i++)
    {
        std::vector<double> x0_vec = {x1_0_diverge[i], x2_0_diverge[i]};
        auto dx = rungeKutta4(f_sistema, x0_vec, ti, tf_diverge, dt, p, p2);
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
    fprintf(gp, "set xrange [-2:2]\n");
    fprintf(gp, "set yrange [-2:2]\n");
    fprintf(gp, "set key off\n");

    std::cout << "Generazione curve di livello...\n";
    // Disegno le curve di livello
    // ==========================================
    // CURVE DI LIVELLO
    //   V(x)     = x1 - 2*(1+ln(x1/2)) + x2 - 2*(1+ln(x2/1))
    // ==========================================
    fprintf(gp, "V(x,y)    = 0.5*(x*x + y*y)\n");

    // Risoluzione della griglia usata per calcolare i contorni (piu' alta = curva piu' liscia)
    fprintf(gp, "set isosamples 300,300\n");
    fprintf(gp, "set contour base\n");
    fprintf(gp, "unset surface\n");

    // --- Curva V(x) ---
    fprintf(gp, "set cntrparam levels discrete 0.25,0.5,2\n");
    fprintf(gp, "set table 'v_level.dat'\n");
    fprintf(gp, "splot [-2:2][-2:2] V(x,y)\n");
    fprintf(gp, "unset table\n");

    // Ripristina lo stato per il grafico 2D finale
    fprintf(gp, "unset contour\n");
    fprintf(gp, "set surface\n");

    // ==========================================
    // COSTRUZIONE DEL COMANDO PLOT E INVIO DATI INLINE
    // ==========================================
    std::cout << "Generazione delle traiettorie...\n";
    fprintf(gp, "plot ");
    for (size_t i = 0; i < tutte_le_traiettorie.size(); ++i)
    {
        // Prepariamo Gnuplot a ricevere x1_0.size() blocchi di dati inline ('-')
        fprintf(gp, "'-' with lines lc rgb 'blue' lw 1.2 notitle, ");
    }
    fprintf(gp, "'v_level.dat' with lines lc rgb 'red' lw 1.5 notitle\n");

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

    if (std::remove("v_level.dat") != 0)
        std::cerr << "Attenzione: impossibile rimuovere v_level.dat" << std::endl;

    std::cout << "Grafico generato con successo!" << std::endl;
    return 0;
}