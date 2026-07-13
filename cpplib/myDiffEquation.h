#include <vector>
#include <array>
#include <complex>
#include <cmath>

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

// Moltiplica matrice * vettore
std::vector<double> matVec(const std::vector<std::vector<double>> &M,
                           const std::vector<double> &v)
{
    size_t rows = M.size(), cols = v.size();
    std::vector<double> result(rows, 0.0);
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            result[i] += M[i][j] * v[j];
    return result;
}

// Somma due vettori
std::vector<double> vecAdd(const std::vector<double> &a,
                           const std::vector<double> &b)
{
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++)
        result[i] = a[i] + b[i];
    return result;
}

// Moltiplica vettore per scalare
std::vector<double> vecScale(const std::vector<double> &v, double s)
{
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); i++)
        result[i] = v[i] * s;
    return result;
}

std::vector<double> lsim(
    const std::vector<std::vector<double>> &A,
    const std::vector<std::vector<double>> &B,
    const std::vector<std::vector<double>> &C,
    const std::vector<std::vector<double>> &D,
    const std::vector<double> &u, // segnale di ingresso
    const std::vector<double> &t) // vettore tempi
{
    size_t n = A.size();     // dimensione stato
    double dt = t[1] - t[0]; // passo temporale

    std::vector<double> x(n, 0.0); // stato iniziale = 0
    std::vector<double> y;         // output

    for (size_t k = 0; k < t.size(); k++)
    {
        std::vector<double> u_k = {u[k]};

        // y(k) = C*x + D*u
        std::vector<double> Cx = matVec(C, x);
        std::vector<double> Du = matVec(D, u_k);
        y.push_back(Cx[0] + Du[0]);

        // x(k+1) = x(k) + dt * (A*x + B*u)   [Eulero]
        std::vector<double> Ax = matVec(A, x);
        std::vector<double> Bu = matVec(B, u_k);
        std::vector<double> dx = vecAdd(Ax, Bu);
        x = vecAdd(x, vecScale(dx, dt));
    }

    return y;
}

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