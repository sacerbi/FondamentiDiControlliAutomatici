#include <vector>
#include <array>
#include <complex>

// Restituisce i due autovalori (eventualmente complessi) di una matrice 2x2
std::array<std::complex<double>, 2> eigenvalues2x2(
    const std::array<std::array<double, 2>, 2> &A)
{
    // Caratteristica: λ² - tr(A)·λ + det(A) = 0
    double tr = A[0][0] + A[1][1];                      // traccia
    double det = A[0][0] * A[1][1] - A[0][1] * A[1][0]; // determinante
    double discriminant = tr * tr - 4.0 * det;

    std::complex<double> sqrt_disc = (discriminant >= 0)
                                         ? std::complex<double>(std::sqrt(discriminant), 0.0)
                                         : std::complex<double>(0.0, std::sqrt(-discriminant));

    return {
        (tr + sqrt_disc) / 2.0,
        (tr - sqrt_disc) / 2.0};
}

// --------------------------
// Operatori utili
// --------------------------
std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] + b[i];
    return r;
}
std::vector<double> operator*(const std::vector<double> &a, double s)
{
    std::vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        r[i] = a[i] * s;
    return r;
}

// --------------------------
// Runge–Kutta 4° ordine
// --------------------------
template <typename Func, typename Params, typename EditableParams>
std::vector<std::vector<double>> rungeKutta4(Func f, std::vector<double> x0,
                                             double t0, double tf, double dt, const Params &p, EditableParams &ep)
{
    int N = static_cast<int>((tf - t0) / dt);
    std::vector<std::vector<double>> result;
    result.reserve(N + 1);

    std::vector<double> x = x0;
    double t = t0;
    result.push_back(x);

    for (int i = 0; i < N; ++i)
    {
        auto k1 = f(t, x, p, ep);
        auto k2 = f(t + dt / 2, x + k1 * (dt / 2), p, ep);
        auto k3 = f(t + dt / 2, x + k2 * (dt / 2), p, ep);
        auto k4 = f(t + dt, x + k3 * dt, p, ep);

        x = x + (k1 + k2 * 2 + k3 * 2 + k4) * (dt / 6);
        t += dt;
        result.push_back(x);
    }

    return result;
}