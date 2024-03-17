#include "hw4.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <string>
#include <utility>
#include "hw2.cpp"
#include <map>

double inf_norm(arma::Mat<double> &P_new, arma::Mat<double> &P_old)
{

    double max_diff = 0;
    for (int i = 0; i < P_new.n_rows; i++)
    {
        for (int j = 0; j < P_new.n_rows; j++)
        {
            double diff = abs(P_new(i, j) - P_old(i, j));
            if (diff > max_diff)
            {
                max_diff = diff;
            }
        }
    }
    return max_diff;
}

void test(string filename)
{

    CNDO2 molecule = CNDO2(filename);
    // print info thint
    const int N = molecule.N;
    arma::Mat<double> P_alpha(N, N, arma::fill::zeros);
    arma::Mat<double> P_beta(N, N, arma::fill::zeros);
    arma::Mat<double> F_alpha(N, N);
    arma::Mat<double> F_beta(N, N);
    arma::Mat<double> C_alpha(N, N);
    arma::Mat<double> C_beta(N, N);

    cout << " Starting Guess P_alpha = P_beta = 0" << endl;
    cout << "P_alpha: " << endl;
    P_alpha.print();
    cout << "P_beta: " << endl;
    P_beta.print();

    int itt = 0;
    double max_mag_change = 1;

    while (max_mag_change > 0.5)
    {

        molecule.calc_F(F_alpha, P_alpha);
        molecule.calc_F(F_beta, P_beta);

        cout << "F_alpha: " << endl;
        F_alpha.print();
        cout << "F_beta: " << endl;
        F_beta.print();

        molecule.calc_C(C_alpha, F_alpha);
        molecule.calc_C(C_beta, F_beta);

        cout << "C_alpha: " << endl;
        C_alpha.print();
        cout << "C_beta: " << endl;
        F_beta.print();

        arma::Mat<double> new_P_alpha(N, N) = P_alpha;
        rma::Mat<double> new_P_beta(N, N) = P_beta;
        molecule.calc_P(new_P_alpha, C_alpha, alpha_);
        molecule.calc_P(new_P_beta, C_beta, beta_);

        cout << "new P_alpha: " << endl;
        new_P_alpha.print();
        cout << "new P_beta: " << endl;
        new_P_beta.print();

        max_mag_change = 0;
    }
}

int main()
{
    test("H2.txt");

    return 0;
}