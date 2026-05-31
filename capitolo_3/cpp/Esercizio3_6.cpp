#include <iostream>
#include <iomanip>
#include <string>
#include "myDiffEquation.h"

std::vector<double> componiPolinomio(const std::vector<int> &i, const std::vector<double> &vMin, const std::vector<double> &vMax)
{
    std::vector<double> v;
    for (int id = 0; id < i.size(); id++)
    {
        int coeff = i[id];
        v.push_back((1 - coeff) * vMin[id] + coeff * vMax[id]);
    }
    return v;
}

void stampaPoli(const std::vector<double> &v, std::string name)
{
    std::cout << "Polinomio " << name << ": [ ";
    for (auto c : v)
        std::cout << c << " ";
    std::cout << "]" << std::endl;
}

int main()
{
    std::vector<double> phi_min = {3.0, 2.0, 6.0, 1.0};
    std::vector<double> phi_max = {4.0, 3.0, 6.0, 2.0};

    int order = phi_max.size();

    std::vector<double> Pa, Pb, Pc, Pd;
    std::vector<int> idx;
    for (int i = 0; i < order; i++)
        idx.push_back(i);

    std::vector<int> idx_a, idx_b, idx_c, idx_d;
    for (auto i : idx)
    {
        idx_a.push_back(i % 4 < 2);
        idx_b.push_back(i % 4 >= 2);
        idx_c.push_back(i % 4 == 0 || i % 4 == 3);
        idx_d.push_back(i % 4 == 1 || i % 4 == 2);
    }

    Pa = componiPolinomio(idx_a, phi_min, phi_max);
    Pb = componiPolinomio(idx_b, phi_min, phi_max);
    Pc = componiPolinomio(idx_c, phi_min, phi_max);
    Pd = componiPolinomio(idx_d, phi_min, phi_max);

    stampaPoli(Pa, "Pa");
    stampaPoli(Pb, "Pb");
    stampaPoli(Pc, "Pc");
    stampaPoli(Pd, "Pd");

    bool aStabile = isStabile(routhTable(Pa));
    bool bStabile = isStabile(routhTable(Pb));
    bool cStabile = isStabile(routhTable(Pc));
    bool dStabile = isStabile(routhTable(Pd));

    std::cout << (aStabile && bStabile && cStabile && dStabile ? "Il polinomio e' stabile" : "Il polinomio NON e' stabile") << std::endl;

    return 0;
}