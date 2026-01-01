#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include "matplotlibcpp.h"
#include "myDiffEquation.h"

namespace plt = matplotlibcpp;

// --------------------------
// Parametri di sistema
// --------------------------
struct Params {
    double Mbar;
    double hbar;
    double kbar;
    double Fe;
    double sbar;
    double mu;
    double eta;
};

// --------------------------
// Equazione differenziale - Controllore ad anello chiuso
// --------------------------
std::vector<double> f_evoluzioneDiscreta(double y0, double w, double mu, int k) {
    std::vector<double> y(k+1);
    y[0] = y0;
    for (int i = 0; i < k; ++i) {
        y[i+1] = mu * (w - y[i]);
    }    
    return y;
}

// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double y0 = 0, k = 4;
    std::vector<double> mus_prop = {-1.1, -0.5, 0.0, 0.5, 1.1};
    std::vector<double> ws_prop = {1.0, 2.0, 3.0};
    std::vector<double> time(k + 1);
    for (int i = 0; i <= k; ++i) time[i] = i;

    // -----------------------------------------
    // Figura 1.7 - Evoluzione del processo, sotto azione del controllore.
    // -----------------------------------------
        for (int w_idx = 0; w_idx < ws_prop.size(); ++w_idx) {
        double w = ws_prop[w_idx];
        std::map<double, std::vector<double>> evoluzioneProcesso;;
        for (auto mu : mus_prop) {
            auto traj = f_evoluzioneDiscreta(y0, w, mu, k);
            evoluzioneProcesso[mu] = traj;
        }
        plt::figure(w_idx + 1);
        for (auto& [mu, C] : evoluzioneProcesso)
            plt::plot(time, C, {{"label", "μ = " + std::to_string(mu)}});
        plt::plot(time, std::vector<double>(time.size(), w), {{"linestyle", "--"}});
        plt::plot(time, std::vector<double>(time.size(), -w), {{"linestyle", "--"}});
        plt::xlabel("Tempo (k)");
        plt::ylabel("Evoluzione y");
        plt::title("Evoluzione del processo per w = " + std::to_string(w));
        plt::legend();
        plt::grid(true);
        plt::show();
    }
    return 0;
}
