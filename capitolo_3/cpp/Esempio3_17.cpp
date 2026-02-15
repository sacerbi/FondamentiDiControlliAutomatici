#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>

// ----------------------------------
// Funzione ausiliaria per verificare che una riga sia nulla
// ----------------------------------
bool isRowZero(const std::vector<double>& v, double epsilon = 0.0001){
    for (auto val : v){
        if(std::abs(val) > epsilon) return false;
    }
    return true;
}

// --------------------------
// Implementazione del sistema differenziale
// --------------------------
std::vector<std::vector<double>> routhTable(const std::vector<double>& poly, double epsilon = 0.0001) {
    int nCoeff = poly.size();
    std::cout << "Routh Table method for " << nCoeff-1 << " degree polynomial. Using epsilon = "<< epsilon << "\n";

    // Creazione della tabella vuota
    std::vector<std::vector<double>> rt(nCoeff, std::vector<double>(static_cast<int>(ceil(nCoeff/2.0)), 0.0));

    // Popolo le prime due righe
    for (int i = 0; i < nCoeff; i++){
        int r = i%2; //Decido se prima o seconda riga a seconda che l'esponente sia pari o dispari
        int c = static_cast<int>(ceil(i/2)); //Procedo con le colonne
        rt[r][c] = poly[i];
    }

    int rowsToDo = nCoeff - 2;
    std::vector<int> maxIndexColumn(rowsToDo, 0);
    for (int i = 0; i < rowsToDo; i++)
        maxIndexColumn[rowsToDo-i-1] = static_cast<int>(ceil(i/2)) + 1;

    for (int i = 2; i < nCoeff; i++){
        if(isRowZero(rt[i-1])){
            std::cout << "\nCaso speciale riga nulla\n";
            int orderAuxPoly = nCoeff - i + 1; //Ordine del polinomio ausiliario
            int nCoeffAux = static_cast<int>(ceil(orderAuxPoly/2.0)); //Num. coeff. polinomio ausiliario
            std::vector<double> auxPoly(nCoeffAux); //Prendo il polinomio ausiliario da usare
            for (int p = 0; p < nCoeffAux; p++){
                auxPoly[p] = rt[i-2][p];
            }
            std::vector<int> auxPow; //Il polinomio ha potenze pari
            for (int k = orderAuxPoly; k >= 0; k-=2){
                auxPow.push_back(k);
            }
            //Assegno la derivata
            for (int k = 0; k < nCoeffAux; k++){
                rt[i-1][k] = auxPoly[k] * auxPow[k];
            }
        }
        else if (std::abs(rt[i-1][0]) < 0.0001){
            std::cout << "\nCaso speciale primo termine nullo\n";
            bool allPositive = true;
            for (int p = 0; p < i-1; p++){
                if (!allPositive) break;
                allPositive &= rt[p][0] >= 0;
            }
            rt[i-1][0] = allPositive ? epsilon : -epsilon;
        }
        for (int j = 0; j < maxIndexColumn[i-2]; j++){
        double h1 = rt[i-2][0];
        double k1 = rt[i-1][0];
        double hi1 = rt[i-2][j+1];
        double ki1 = rt[i-1][j+1];
        rt[i][j] = hi1 - ((h1*ki1)/k1);
    }   
    }
     
    return rt;
}

bool isStabile(const std::vector<std::vector<double>>& rt) {
    bool tuttiPositivi = true;
    bool tuttiNegativi = true;
    for (const auto& riga : rt) {
        if (riga[0] <= 0) tuttiPositivi = false;
        if (riga[0] >= 0) tuttiNegativi = false;
    }
    return tuttiPositivi || tuttiNegativi;
}

std::vector<double> componiPolinomio(const std::vector<int>& i, const std::vector<double>& vMin, const std::vector<double>& vMax){
    std::vector<double> v;
    for (int id = 0; id < i.size(); id++){
        int coeff = i[id];
        v.push_back((1-coeff)*vMin[id] + coeff * vMax[id]);
    }
    return v;
}

void stampaPoli(const std::vector<double>& v, std::string name){
    std::cout << "Polinomio " << name << ": [ ";
    for(auto c : v) std::cout << c << " ";
    std::cout << "]" << std::endl;
}


int main() {
    std::vector<double> phi_min = {0.9, 13.5, 76.5, 202.5, 246.6, 108};
    std::vector<double> phi_max = {1.1, 16.5, 93.5, 247.5, 301.4, 132};

    int order = phi_max.size();

    std::vector<double> Pa, Pb, Pc, Pd;
    std::vector<int> idx;
    for (int i = 0; i < order; i++) idx.push_back(i);

    std::vector<int> idx_a, idx_b, idx_c, idx_d;
    for(auto i : idx){
        idx_a.push_back(i%4 < 2);
        idx_b.push_back(i%4 >= 2);
        idx_c.push_back(i%4 == 0 || i%4 == 3);
        idx_d.push_back(i%4 == 1 || i%4 == 2);
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

    std::cout << (aStabile && bStabile && cStabile && dStabile ? "Il polinomio è stabile" : "Il polinomio NON è stabile") << std::endl;

    return 0;
}