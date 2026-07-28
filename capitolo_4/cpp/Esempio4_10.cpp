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
struct EmptyParams
{
};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_sistema(double t, const std::vector<double> &x, const EmptyParams &p, const EmptyParams &p2)
{
    std::vector<double> dx(2);
    dx[0] = -x[1];
    dx[1] = x[0] * x[0] * x[1] + x[0] - x[1];
    return dx;
}

// --------------------------
// MAIN
// --------------------------
int main()
{
    // Parametri di simulazione
    EmptyParams p;
    double dt = 0.1;
    double ti = 0.0, tf = 20.0, tf_diverge = 2.2;

    // Vettori degli stati iniziali
    std::vector<double> x1_0 = {-2.0, 2.0, -2.0, 2.0};
    std::vector<double> x2_0 = {0.0, 0.0, 0.1, -0.1};
    std::vector<double> x1_0_diverge = {-2.0, 2.0};
    std::vector<double> x2_0_diverge = {-1.0, 1.0};

    // Struttura per memorizzare tutte le traiettorie in RAM prima di inviarle a Gnuplot
    std::vector<std::vector<std::vector<double>>> tutte_le_traiettorie;

    std::cout << "Calcolo delle traiettorie in corso...\n";
    // Simulazione e salvataggio in memoria
    for (size_t i = 0; i < x1_0.size(); i++)
    {
        std::vector<double> x0_vec = {x1_0[i], x2_0[i]};
        auto dx = rungeKutta4(f_sistema, x0_vec, ti, tf, dt, p, p);
        tutte_le_traiettorie.push_back(dx);
    }
    for (size_t i = 0; i < x1_0_diverge.size(); i++)
    {
        std::vector<double> x0_vec = {x1_0_diverge[i], x2_0_diverge[i]};
        auto dx = rungeKutta4(f_sistema, x0_vec, ti, tf_diverge, dt, p, p);
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
    fprintf(gp, "set xrange [-4:4]\n");
    fprintf(gp, "set yrange [-4:4]\n");
    fprintf(gp, "set key off\n");

    // ==========================================
    // CURVE DI LIVELLO
    //   V(x)     = 0.5*(3*x1^2 - 2*x1*x2 + 2*x2^2)   -> livello 2.2
    //   Vdot(x)  = -x1^2 -x2^2 -x1^3*x2 + 2*x1^2*x2^2 -> livello 0
    // ==========================================
    fprintf(gp, "V(x,y)    = 0.5*(3*x*x - 2*x*y + 2*y*y)\n");
    fprintf(gp, "Vdot(x,y) = -x*x - y*y - x*x*x*y + 2*x*x*y*y\n");

    // Risoluzione della griglia usata per calcolare i contorni (piu' alta = curva piu' liscia)
    fprintf(gp, "set isosamples 300,300\n");
    fprintf(gp, "set contour base\n");
    fprintf(gp, "unset surface\n");

    // --- Curva V(x) = 2.2 ---
    fprintf(gp, "set cntrparam levels discrete 2.2\n");
    fprintf(gp, "set table 'v_level.dat'\n");
    fprintf(gp, "splot [-4:4][-4:4] V(x,y)\n");
    fprintf(gp, "unset table\n");

    // --- Curva Vdot(x) = 0 ---
    fprintf(gp, "set cntrparam levels discrete 0\n");
    fprintf(gp, "set table 'vdot_level.dat'\n");
    fprintf(gp, "splot [-4:4][-4:4] Vdot(x,y)\n");
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
        // Prepariamo Gnuplot a ricevere un blocco di dati inline ('-') per ogni traiettoria
        fprintf(gp, "'-' with lines lc rgb 'blue' lw 1.2 notitle, ");
    }
    fprintf(gp, "'v_level.dat' with lines lc rgb 'red' lw 1.5 notitle, ");
    fprintf(gp, "'vdot_level.dat' with lines lc rgb 'green' dashtype 2 lw 1.5 notitle\n");

    // Invio sequenziale dei blocchi di dati delle traiettorie
    // (l'ordine deve rispettare quello dei placeholder '-' nel comando plot sopra)
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
    if (std::remove("vdot_level.dat") != 0)
        std::cerr << "Attenzione: impossibile rimuovere vdot_level.dat" << std::endl;

    std::cout << "Grafico generato con successo!" << std::endl;
    return 0;
}