#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <string>
#include <utility>

// namespace hw3
//{
using namespace std;

using namespace arma;

//}

// using namespace hw3;

struct Function
{
    vector<double> R;
    vector<double> norms;
    vector<double> coefs;
    vector<double> exponets;
    int l;
    int m;
    int n;
    double int_pot;
    vector<int> orbital;

    Function(vector<double> R, vector<double> norms, vector<double> coefs, vector<double> exponets, int l, int m, int n, double int_pot, vector<int> orbital) : R(R), norms(norms), coefs(coefs), exponets(exponets), l(l), m(m), n(n), int_pot(int_pot), orbital(orbital)
    {

        // R.resize(3, 0.0);
        // coefs.resize(3, 0.0);
        // norms.resize(3, 0.0);
        // l = 0;
        // m = 0;
        // n = 0;
        // int_pot = 0.0;
    }
};

class huckle
{

public:
    huckle(string xyz_filepath);
    void read_xyz(string filepath);
    void read_basis(string filepath, int atomic_num);
    // void setup(string xyz_filepath, vector<pair<string, int>> basis_filepaths);
    void calc_basis();
    void calc_S();
    void calc_H();
    void calc_X();
    void calc_HAM();
    void calc_C();
    void calc_e();
    void print_S();
    void print_HAM();
    void print_X();
    void print_H();
    Mat<double> get_S();
    Mat<double> get_H();
    Mat<double> get_X();
    Mat<double> get_HAM();
    Mat<double> get_C();
    Mat<double> get_e();
    double get_total_energy();
    void add_norms();
    void test();

protected:
    int num_atoms;
    vector<double> exponets_H;
    vector<double> exponets_C;
    vector<double> coefs_1s;
    vector<double> coefs_2s;
    vector<double> coefs_3p;
    int charge;
    vector<vector<double>> H_atoms;
    vector<vector<double>> C_atoms;
    int a;
    int b;
    int N;
    int n;
    vector<Function> basis;

    double overlap_integral_1d_analytical(double X_A, double X_B, int l_A, int l_B, double alpha, double beta);
    int factorial(int num);
    int double_factorial(int num);
    double X_p(double X_A, double X_B, double alpha, double beta);
    int binomial_coefficient(int m, int n);

    Mat<double> S;
    Mat<double> H;
    Mat<double> X;
    Mat<double> HAM;
    Mat<double> C;
    Mat<double> V;
    vec e;
};

huckle::huckle(string xyz_filepath)
{
    read_xyz(xyz_filepath);
    a = C_atoms.size();
    b = H_atoms.size();
    N = 4 * a + b;

    try
    {

        if (b % 2 != 0)
        {
            throw runtime_error("Valance electrons must be in pairs");
        }

        // Code here would not execute if an exception is thrown above
    }
    catch (const runtime_error &e)
    {
        cerr << "An error occurred: " << e.what() << endl;
        cerr << " 2a + b/2 for C_aH_b must be an interger value." << endl;
        // The program continues here after handling the exception
    }

    n = 2 * a + (b / 2);
    //}

    X.set_size(N, N);
    S.set_size(N, N);
    H.set_size(N, N);
    HAM.set_size(N, N);
    C.set_size(N, N);

    if (a > 0)
    {
        read_basis("C_STO3G.txt", 6);
    }

    if (b > 0)
    {
        read_basis("H_STO3G.txt", 1);
    }
}
int huckle::factorial(int num)
{
    int product = 1;
    if (num == -1)
    {
    }
    else
    {
        for (int i = 1; i < num + 1; i++)
        {
            product *= i;
        }
    }

    return product;
}

int huckle::double_factorial(int num)
{
    int product = 1;
    if (num % 2 == 0)
    {
        for (int i = 2; i < num + 1; i += 2)
        {
            product *= i;
        }
    }
    else
    {
        for (int i = 1; i < num + 1; i += 2)
        {
            product *= i;
        }
    }
    return product;
}
double huckle::X_p(double X_A, double X_B, double alpha, double beta)
{
    return (alpha * X_A + beta * X_B) / (alpha + beta);
}

int huckle::binomial_coefficient(int m, int n)
{
    return factorial(m) / (factorial(n) * factorial(m - n));
}

double huckle::overlap_integral_1d_analytical(double X_A, double X_B, int l_A, int l_B, double alpha, double beta)
{
    double sum = 0;
    double X_P = (alpha * X_A + beta * X_B) / (alpha + beta);
    for (int i = 0; i < l_A + 1; i++)
    {
        for (int j = 0; j < l_B + 1; j++)
        {
            if ((i + j) % 2 == 0)
            {
                sum += binomial_coefficient(l_A, i) * binomial_coefficient(l_B, j) * double_factorial(i + j - 1) * pow(X_P - X_A, l_A - i) * pow(X_P - X_B, l_B - j) / pow(2 * (alpha + beta), (i + j) / 2);
            }
            else
            {
                continue;
            }
        }
    }
    double result = exp(-alpha * beta * pow(X_A - X_B, 2) / (alpha + beta)) * pow(M_PI / (alpha + beta), 0.5) * sum;
    return result;
}

void huckle::read_xyz(string filepath)
{

    ifstream myfile(filepath);

    if (myfile.is_open())
    {
        myfile >> num_atoms;
        myfile >> charge;

        int atomic_num;
        vector<double> coords(3);

        for (int i = 0; i < num_atoms; i++)
        {

            myfile >> atomic_num;
            myfile >> coords[0];
            myfile >> coords[1];
            myfile >> coords[2];
            if (atomic_num == 1)
            {

                H_atoms.push_back(coords);
            }
            else if (atomic_num == 6)
            {
                C_atoms.push_back(coords);
            }
        }
    }

    myfile.close();
}

void huckle::read_basis(string filepath, int atomic_num)
{
    fstream myfile(filepath);

    if (myfile.is_open())
    {
        double exponet;
        double coef;
        if (atomic_num == 1)
        {
            double coef;

            for (int i = 0; i < 3; i++)
            {
                myfile >> exponet;
                myfile >> coef;
                exponets_H.push_back(exponet);
                coefs_1s.push_back(coef);
            }
        }
        if (atomic_num == 6)
        {

            for (int i = 0; i < 3; i++)
            {
                myfile >> exponet;
                exponets_C.push_back(exponet);
                myfile >> coef;
                coefs_2s.push_back(coef);
                myfile >> coef;
                coefs_3p.push_back(coef);
            }
        }
    }
}

void huckle::calc_basis()
{

    for (int i = 0; i < b; i++)
    {

        vector<double> norms;
        vector<int> orb = {0, 0, 0};
        Function fun = Function(H_atoms[i], norms, coefs_1s, exponets_H, 0, 0, 1, -13.6, orb);
        basis.push_back(fun);
    }

    for (int i = 0; i < a; i++)
    {

        vector<double>
            norms;
        vector<int> orb = {0, 0, 0};
        // Function fun = Function(H_atoms[i], norms, coefs_2s, exponets_C, 0, 0, 0, -21.4, orb);
        basis.push_back(fun);

        int count = 0;
        for (int m = -1; m < 2; m++)
        {
            vector<double> norms2;
            vector<int> orb2;

            if (count == 0)
            {
                orb2 = {1, 0, 0};
            }
            else if (count == 1)
            {
                orb2 = {0, 1, 0};
            }

            else if (count == 2)
            {
                orb2 = {0, 0, 1};
            }

            Function fun2 = Function(C_atoms[i], norms2, coefs_3p, exponets_C, 1, 3, m, -11.4, orb2);

            basis.push_back(fun2);

            count++;
        }
    }
}

void huckle::add_norms()
{

    for (int k = 0; k < basis.size(); k++)
    {
        Function &fun = basis[k];
        for (int i = 0; i < 3; i++)
        {
            double overlap = 1;

            for (int j = 0; j < 3; j++)
            {
                overlap *= overlap_integral_1d_analytical(fun.R[j], fun.R[j], fun.orbital[j], fun.orbital[j], fun.exponets[i], fun.exponets[i]);
            }
            fun.norms.push_back(1 / pow(overlap, 0.5));
        }
    }
}

void huckle::calc_S()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {

            S.at(i, j) = 0.0;
            Function &u = basis[i];
            Function &v = basis[j];

            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    double overlap = 1;

                    for (int z = 0; z < 3; z++)
                    {
                        overlap *= overlap_integral_1d_analytical(u.R[z], v.R[z], u.orbital[z], v.orbital[z], u.exponets[k], v.exponets[l]);
                    }

                    S.at(i, j) += overlap * u.norms[k] * v.norms[l] * u.coefs[k] * v.coefs[l];
                }
            }
        }
    }
}

void huckle::print_S()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << S.at(i, j) << " ";
        }
        cout << endl;
    }
}

Mat<double> huckle::get_S()
{

    calc_S();
    return S;
}

void huckle::calc_X()
{

    mat P;
    arma::Mat<double> D = arma::zeros(N, N);
    vec v;

    eig_sym(v, P, S);

    for (int i = 0; i < v.size(); i++)
    {
        if (abs(v(i)) > numeric_limits<double>::epsilon())
        {
            D.at(i, i) = 1 / pow(abs(v(i)), 0.5);
        }
    }

    X = P * D * P.t();
}

void huckle::print_X()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << X(i, j) << " ";
        }
        cout << endl;
    }
}

Mat<double> huckle::get_X()
{
    calc_X();
    return X;
}

void huckle::calc_H()
{

    for (int i = 0; i < N; i++)
    {
        Function &u = basis[i];
        H.at(i, i) = u.int_pot;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)

        {
            Function &k = basis[i];
            Function &v = basis[j];
            if (j != i)
            {

                H.at(i, j) = (1.75 / 2) * (H(i, i) + H(j, j)) * S(i, j);
            }
            else
            {
                continue;
            }
        }
    }
}

Mat<double> huckle::get_H()
{
    calc_H();
    return H;
}

void huckle::print_H()
{
    // for (int i = 0; i < N; i++)
    // {
    //     for (int j = 0; j < N; j++)
    //     {
    //         cout << H.at(i, j) << " ";
    //     }
    //     cout << endl;
    // }
    H.print();
}

void huckle::calc_HAM()
{
    HAM = X * H * X.t();
}

Mat<double> huckle::get_HAM()
{
    calc_HAM();
    return HAM;
}

void huckle::calc_e()
{

    arma::eig_sym(e, V, HAM);
}

void huckle::calc_C()
{
    C = X * V;
}

void huckle::print_HAM()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout << HAM(i, j) << " ";
        }
        cout << endl;
    }
}

Mat<double> huckle::get_C()
{
    calc_C();
    return C;
}

double huckle::get_total_energy()
{
    double E = 0.0;
    for (int i = 0; i < n; i++)
    {
        E += 2 * e(i);
    }
    return E;
}

// void huckle::test()
// {

//     arma::mat22 ham;
//     ham << -13.6000 << -15.7050 << arma::endr
//         << -15.7050 << -13.6000 << arma::endr;

//     mat h;
//     h = X * ham * X.t();
//     h.print();
// }

// int main()
// {
//     // #1

//     huckle H2 = huckle("H2.txt");
//     H2.calc_basis();
//     H2.add_norms();
//     H2.calc_S();
//     // H2.print_S();
//     H2.calc_X();
//     H2.calc_H();
//     cout << " Hamiltonian Matrix for H2" << endl;
//     H2.print_H();
//     // H2.print_X();
//     H2.calc_HAM();
//     // H2.print_HAM();

//     H2.calc_e();

//     // H2.test();
//     // H2.print_H();

//     double H2_energy = H2.get_total_energy();
//     cout << "The total energy of H2 is: " << H2_energy << " eV." << endl;

//     double H2_bond_energy = H2_energy - (-27.2);
//     cout << "The bond energy of H2 is: " << H2_bond_energy << " eV." << endl;

//     // #2
//     huckle C2H2 = huckle("C2H2.txt");
//     C2H2.calc_basis();
//     C2H2.add_norms();
//     C2H2.calc_S();
//     C2H2.calc_X();
//     C2H2.calc_H();
//     cout << " Hamiltonian Matrix for C2H2" << endl;
//     C2H2.print_H();
//     C2H2.calc_HAM();
//     C2H2.calc_e();

//     huckle C2H4 = huckle("C2H4.txt");
//     C2H4.calc_basis();
//     C2H4.add_norms();
//     C2H4.calc_S();
//     C2H4.calc_X();
//     C2H4.calc_H();
//     cout << " Hamiltonian Matrix for C2H4" << endl;
//     C2H4.print_H();
//     C2H4.calc_HAM();
//     C2H4.calc_e();

//     cout << endl;
//     double C2H2_energy = C2H2.get_total_energy();
//     cout << "The total energy of C2H2 is: " << C2H2_energy << " eV." << endl;

//     double C2H4_energy = C2H4.get_total_energy();
//     cout << "The total energy of C2H4 is: " << C2H4_energy << " eV." << endl;

//     double change_e = (C2H4_energy - C2H2_energy - H2_energy) * 96.485;
//     cout << "The total change in ehtralopy for the reaction C2H2 + H2 = C2H4 is: " << C2H4_energy << " kJ/mol." << endl;

//     return 0;
// }