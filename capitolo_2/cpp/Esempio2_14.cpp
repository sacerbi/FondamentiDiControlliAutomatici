#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <map>
#include "matplotlibcpp.h"
#include "myDiffEquation.h"

namespace plt = matplotlibcpp;
// --------------------------
// MAIN
// --------------------------
int main() {
    // Dati
    double u1 = 1.0/8.0, u2 = 1.0/4.0, u3 = 1.0/2.0;
    double dx = 0.01;
    double x0 = -2, xf = 2;
    int N = static_cast<int>((xf - x0) / dx);
    std::vector<double> x(N + 1);
    for (int i = 0; i <= N; ++i) x[i] = x0 + i * dx;
    // Calcolo il comportamento dei punti di equilibrio al variare di x
    std::vector<double> y1(N+1);
    std::vector<double> y2(N+1);
    std::vector<double> y3(N+1);

    for (int i = 0; i<= N; i++){
        y1[i] = x[i] * x[i] + x[i] + u1;
        y2[i] = x[i] * x[i] + x[i] + u2;
        y3[i] = x[i] * x[i] + x[i] + u3;
    }
    // Plot Figura 2.4
    plt::figure(1);
    plt::plot(x, y1, {{"label", "y(x), ubar = 0.125"}});
    plt::plot(x, y2, {{"label", "y(x), ubar = 0.25"}});
    plt::plot(x, y3, {{"label", "y(x), ubar = 0.5"}});
    plt::plot({-2, 2}, {0, 0}, {{"label", "y = 0"}, {"linestyle", "--"}});
    plt::xlim(-2, 2);
    plt::ylim(-1, 2);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::title("Variazione equilibrio al variare di ubar");
    plt::legend();
    plt::grid(true);
    
    plt::show();

    plt::detail::_interpreter::kill();

    return 0;
}
