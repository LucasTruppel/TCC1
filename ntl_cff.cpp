#include <iostream>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <vector>
using namespace std;
using namespace NTL;

bool DEBUG = false;

ZZ_pE* createAllElementsOfGF(unsigned p, unsigned n) {
    // define GF(p)
    ZZ_p::init(ZZ(p));

    // generate an irreducible polynomial P of degree n over GF(p)
    ZZ_pX P;
    BuildIrred(P, n);

    // define GF(p^n)
    ZZ_pE::init(P);

    if (DEBUG) std::cout << "\nGL elements:"  << std::endl;
    auto* elements = new ZZ_pE[power_long(p, n)];

    // Iterate over all elements in GF(p^n)
    long num_elements = power_long(p, n);  // Total number of elements in GF(p^n)
    for (long i = 0; i < num_elements; ++i) {
        ZZ_pE element;
        ZZ_pX poly;
        long temp = i;
        // Create the polynomial representation of the element
        for (unsigned j = 0; j < n; ++j) {
            SetCoeff(poly, j, temp % p);
            temp /= p;
        }
        element = to_ZZ_pE(poly);
        if (DEBUG) std::cout << "Element: " << element << std::endl;
        elements[i] = element;
    }
    return elements;
}

ZZ_pEX* getPolynomials(unsigned int q, unsigned int k, const ZZ_pE *elements, unsigned int n) {
    auto* polynomials = new ZZ_pEX[n];
    unsigned poly_index = 0;
    std::vector<unsigned> current(k + 1, 0);

    while (true) {
        ZZ_pEX f;
        for (int i = 0; i <= k; ++i) {
            SetCoeff(f, k-i, elements[current[i]]);
        }
        polynomials[poly_index] = f;
        poly_index++;

        // Increment coefficients starting from the lowest degree
        int pos = k;
        while (pos >= 0 && ++current[pos] >= q) {
            current[pos] = 0;
            pos--;
        }
        if (pos < 0) break;
    }
    return polynomials;
}

vector<vector<int>> buildCFF(unsigned q, unsigned k, ZZ_pE* elements) {
    unsigned n = power_long(q, k+1);
    ZZ_pEX* polynomials = getPolynomials(q, k, elements, n);

    vector<vector<int>> matrix(q * q, vector<int>(n, 0));
    int x = 0;
    int fx = 0;
    if (DEBUG) std::cout << "\n CFF matrix creation:" << std::endl;
    for (int i = 0; i < q * q; i++) {
        for (int j = 0; j < n; j++) {
            if (DEBUG) {
                std::cout << "x: " << x << std::endl
                          << "fx: " << fx << std::endl
                          << "polynomials[j]: " << polynomials[j] << std::endl
                          << "elements[x]: " << elements[x] << std::endl
                          << "eval: " << eval(polynomials[j], elements[x]) << std::endl
                          << "elements[fx]: " << elements[fx] << std::endl << std::endl;
            }
            if (eval(polynomials[j], elements[x]) == elements[fx]) {
                matrix[i][j] = 1;
            }
        }
        fx++;
        if (fx == q) {
            fx = 0;
            x++;
        }
    }

    delete[] polynomials;
    return matrix;
}

void printMatrix(unsigned int numLines, unsigned int numColumns, const vector<vector<int>> &matrix) {
    for (int i = 0; i < numLines; i++) {
        for (int j = 0; j < numColumns; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    unsigned p = 3;
    unsigned n = 1;
    unsigned k = 2;
    unsigned q = power_long(p, n);
    std::cout << "p: " << p <<  " n: " << n << " k: " << k << std::endl;

    ZZ_pE* elements = createAllElementsOfGF(p, n);
    vector<vector<int>> CFF = buildCFF(power_long(p, n), k, elements);
    std::cout << "CFF: " << std::endl;
    printMatrix(q * q, power_long(q, k+1), CFF);

    delete[] elements;
    return 0;
}
