#include <iostream>
#include <iomanip>
#include <cmath>
#include "myDiffEquation.h"
#include <fstream> // <--- AGGIUNTO per salvare i file
#include <cstdlib> // <--- AGGIUNTO per lanciare Gnuplot

// --------------------------
// Parametri del sistema
// --------------------------
struct EmptyParams{};

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<double> f_sistema(double t, const std::vector<double> &x, const EmptyParams &p, const EmptyParams &p2)
{
    std::vector<double> dx(2);
    dx[0] = -std::pow(x[0],3) + x[0]*std::pow(x[1],2) - std::pow(x[1],2);
    dx[1] = x[0]*x[1] - std::pow(x[0],2)*x[1] - std::pow(x[1],3);
    return dx;
}

int main()
{
    // Parametri di simulazione
    EmptyParams p;
    double dt = 0.01;
    double ti = 0.0, tf = 25.0;
    std::vector<double> x0 = {2.0, 1.0};
    int N = static_cast<int>((tf - ti) / dt);
    std::vector<double> time(N + 1);
    for (int i = 0; i <= N; ++i)
        time[i] = ti + i * dt;

    // Vettori degli stati iniziali (presi dal grafico riportato a pag. 103)
    std::vector<double> x1_0 = {-3.0, -3.0, -3.0, -3.0, -3.0, -1.5, -1.5, -1.0, -1.0,  0.0, 0.0,  1.0, 1.0,  2.0, 2.0,  3.0,  3.0, 3.0, 3.0, 3.0};
    std::vector<double> x2_0 = {-3.0, -1.0, 0.0, 1.0, 3.0, -3.0, 3.0, -3.0, 3.0,  -3.0, 3.0,  -3.0, 3.0,  -3.0, 3.0,  -3.0,  -1.0, 0.0, 1.0, 3.0};

    // Apre file per salvare le traiettorie (verrà letto da Gnuplot)
    std::string filename = "esempio4_6.dat";
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Errore nell'apertura del file dati per Gnuplot." << std::endl;
        return 1;
    }

    std::cout << "Calcolo delle traiettorie in corso...\n";
    // Simulazione e calcolo delle traiettorie
    for(size_t i = 0; i < x1_0.size(); i++){
        auto x0 = {x1_0[i], x2_0[i]};
        auto dx = rungeKutta4(f_sistema, x0, ti, tf, dt, p, p);

        // Scrittura nel file formato "blocchi Gnuplot" 
        // (ogni blocco è separato da due righe vuote)
        for (const auto &x : dx) {
            outfile << x[0] << " " << x[1] << "\n";
        }
        outfile << "\n\n";
    }
    outfile.close();
    std::cout << "Dati generati con successo. Avvio di Gnuplot..." << std::endl;

    

    // 4. Avvio ed esecuzione di Gnuplot via Pipe
    // L'opzione -persistent lascia aperta la finestra del grafico al termine del programma C++
    FILE *gp = popen("gnuplot -persistent", "w");
    if (!gp) {
        std::cerr << "Errore: impossibile avviare Gnuplot. Assicurati che sia installato e nel PATH di sistema." << std::endl;
        return 1;
    }

    // Impostazioni estetiche del grafico
    fprintf(gp, "set title 'Traiettorie e linee di livello (di Lyapunov)' font ',12'\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set size square\n"); // Mantiene le proporzioni 1:1 come 'axis square' di Matlab
    fprintf(gp, "set xrange [-3:3]\n");
    fprintf(gp, "set yrange [-3:3]\n");
    fprintf(gp, "set key off\n");     // Nasconde la legenda per non pasticciare il grafico

    // Disegno le curve di livello di V = 0.5*(x^2 + y^2) come circonferenze rosse
    // Se V = c  ->  x^2 + y^2 = 2*c  ->  raggio = sqrt(2*c)
    std::vector<double> livelli_v = {0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
    for (size_t i = 0; i < livelli_v.size(); ++i) {
        double raggio = std::sqrt(2.0 * livelli_v[i]);
        fprintf(gp, "set object %zu circle at 0,0 size %f fc rgb 'red' fs empty lw 1.5\n", i + 1, raggio);
    }

    // Costruzione del comando per plottare tutti i 20 blocchi dal file
    fprintf(gp, "plot ");
    for (size_t i = 0; i < x1_0.size(); ++i) {
        // 'index i' indica a gnuplot di leggere solo l'i-esimo blocco di dati
        fprintf(gp, "'%s' index %zu with lines lc rgb 'blue' lw 1.2", filename.c_str(), i);
        if (i < x1_0.size() - 1) {
            fprintf(gp, ", ");
        }
    }
    fprintf(gp, "\n");

    fflush(gp); // Forza l'invio dei comandi al processo Gnuplot
    pclose(gp); // Chiude la pipe

    std::cout << "Grafico generato con successo!" << std::endl;

    std::remove("esempio4_6.dat");
    return 0;
}