#include <iostream>
#include <cstdio>
#include <cmath>

// Gestione della compatibilità cross-platform per le pipe
#ifdef _WIN32
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

int main() {
    // Apre una pipe verso Gnuplot (-persist mantiene le finestre aperte alla fine del programma)
    FILE* gp = POPEN("gnuplot -persist", "w");
    if (!gp) {
        std::cerr << "Errore: Impossibile trovare Gnuplot. Assicurati che sia installato e nel PATH." << std::endl;
        return 1;
    }

    // ==========================================
    // FIGURA 1: Z1 = 3*X^2 + 5*Y^2
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura 1'\n"); // Apre la prima finestra
    fprintf(gp, "set xrange [-2:2]\n");
    fprintf(gp, "set yrange [-2:2]\n");
    fprintf(gp, "set zrange [0:35]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set zlabel 'W'\n");
    fprintf(gp, "set view 60, 322\n"); // Corrisponde approssimativamente a view(3) di Octave
    fprintf(gp, "splot '-' with lines title 'Z1 = 3*x^2 + 5*y^2'\n");

    // Simulazione di meshgrid e calcolo della superficie
    for (double x = -2.0; x <= 2.01; x += 0.2) {
        for (double y = -2.0; y <= 2.01; y += 0.2) {
            double z = 3.0 * (x * x) + 5.0 * (y * y);
            fprintf(gp, "%f %f %f\n", x, y, z);
        }
        fprintf(gp, "\n"); // RIGA VUOTA: Fondamentale per far capire a Gnuplot la struttura a griglia (mesh)
    }
    fprintf(gp, "e\n"); // 'e' indica a Gnuplot la fine dei dati per questo grafico


    // ==========================================
    // FIGURA 2: Z = X^2 + Y^2 - X^2*Y
    // ==========================================
    fprintf(gp, "set terminal qt 2 title 'Figura 2'\n"); // Apre la seconda finestra
    fprintf(gp, "set xrange [-2:2]\n");
    fprintf(gp, "set yrange [-2:2]\n");
    fprintf(gp, "set zrange [0:10]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set zlabel 'W'\n");
    fprintf(gp, "set view 60, 322\n");
    fprintf(gp, "splot '-' with lines title 'Z = x^2 + y^2 - x^2*y'\n");

    for (double x = -2.0; x <= 2.01; x += 0.2) {
        for (double y = -2.0; y <= 2.01; y += 0.2) {
            double z = (x * x) + (y * y) - (x * x * y);
            fprintf(gp, "%f %f %f\n", x, y, z);
        }
        fprintf(gp, "\n"); 
    }
    fprintf(gp, "e\n");


    // ==========================================
    // FIGURA 3: Z = X^2
    // ==========================================
    fprintf(gp, "set terminal qt 3 title 'Figura 3'\n"); // Apre la terza finestra
    fprintf(gp, "set xrange [-2:2]\n");
    fprintf(gp, "set yrange [-2:2]\n");
    fprintf(gp, "set zrange [0:4]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set zlabel 'W'\n");
    fprintf(gp, "set view 60, 322\n");
    fprintf(gp, "splot '-' with lines title 'Z = x^2'\n");

    for (double x = -2.0; x <= 2.01; x += 0.2) {
        for (double y = -2.0; y <= 2.01; y += 0.2) {
            double z = 1.0 * (x * x);
            fprintf(gp, "%f %f %f\n", x, y, z);
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");

    // Chiude la pipe
    PCLOSE(gp);
    return 0;
}