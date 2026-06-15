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
    // FIGURA 1: Z1 = alpha*(1 - cos X) + beta*Y^2
    // ==========================================
    fprintf(gp, "set terminal qt 1 title 'Figura 1'\n"); // Apre la prima finestra
    fprintf(gp, "set xrange [-4:4]\n");
    fprintf(gp, "set yrange [-4:4]\n");
    fprintf(gp, "set zrange [0:6]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set zlabel 'W'\n");
    fprintf(gp, "set view 60, 322\n"); // Corrisponde approssimativamente a view(3) di Octave
    fprintf(gp, "splot '-' with lines title 'Z1 = (1-cos x_1) + x_2^2'\n");

    double alpha = 1, beta = 1;
    // Simulazione di meshgrid e calcolo della superficie
    for (double x = -4.0; x <= 4.01; x += 0.2) {
        for (double y = -4.0; y <= 4.01; y += 0.2) {
            double z = alpha * (1 - std::cos(x)) + beta * (y * y);
            fprintf(gp, "%f %f %f\n", x, y, z);
        }
        fprintf(gp, "\n"); // RIGA VUOTA: Fondamentale per far capire a Gnuplot la struttura a griglia (mesh)
    }
    fprintf(gp, "e\n"); // 'e' indica a Gnuplot la fine dei dati per questo grafico


    // ==========================================
    // FIGURA 2: Z = e^(X^2+Y^2) - 1
    // ==========================================
    fprintf(gp, "set terminal qt 2 title 'Figura 2'\n"); // Apre la seconda finestra
    fprintf(gp, "set xrange [-1:1]\n");
    fprintf(gp, "set yrange [-1:1]\n");
    fprintf(gp, "set zrange [0:7]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set zlabel 'W'\n");
    fprintf(gp, "set view 60, 322\n");
    fprintf(gp, "splot '-' with lines title 'Z = e^(x_1^2 + x_2^2) -1'\n");

    for (double x = -1.0; x <= 1.01; x += 0.1) {
        for (double y = -1.0; y <= 1.01; y += 0.1) {
            double z = std::exp(x*x + y*y) - 1;
            fprintf(gp, "%f %f %f\n", x, y, z);
        }
        fprintf(gp, "\n"); 
    }
    fprintf(gp, "e\n");

    // ==========================================
    // FIGURA 3: Z = sum(|x_i|)
    // ==========================================
    fprintf(gp, "set terminal qt 3 title 'Figura 3'\n"); // Apre la terza finestra
    fprintf(gp, "set xrange [-1:1]\n");
    fprintf(gp, "set yrange [-1:1]\n");
    fprintf(gp, "set zrange [0:2]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set zlabel 'W'\n");
    fprintf(gp, "set view 60, 322\n");
    fprintf(gp, "splot '-' with lines title 'Z = E(|x_i|)'\n");

    for (double x = -1.0; x <= 1.01; x += 0.1) {
        for (double y = -1.0; y <= 1.01; y += 0.1) {
            double z = std::abs(x) + std::abs(y);
            fprintf(gp, "%f %f %f\n", x, y, z);
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");

    // ==========================================
    // FIGURA 4: Z = max(|x_i|)
    // ==========================================
    fprintf(gp, "set terminal qt 4 title 'Figura 4'\n"); // Apre la quarta finestra
    fprintf(gp, "set xrange [-1:1]\n");
    fprintf(gp, "set yrange [-1:1]\n");
    fprintf(gp, "set zrange [0:1]\n");
    fprintf(gp, "set grid\n");
    fprintf(gp, "set xlabel 'x_1'\n");
    fprintf(gp, "set ylabel 'x_2'\n");
    fprintf(gp, "set zlabel 'W'\n");
    fprintf(gp, "set view 60, 322\n");
    fprintf(gp, "splot '-' with lines title 'Z = max(|x_i|)'\n");

    for (double x = -1.0; x <= 1.01; x += 0.1) {
        for (double y = -1.0; y <= 1.01; y += 0.1) {
            double z = std::max(std::abs(x), std::abs(y));
            fprintf(gp, "%f %f %f\n", x, y, z);
        }
        fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");

    // Chiude la pipe
    PCLOSE(gp);
    return 0;
}