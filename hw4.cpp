#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <string>
#include <utility>
#include <map>
using namespace std;
// using namespace arma;
enum electrons
{
    alpha_,
    beta_
};

struct Function
{
    int atomic_num;
    vector<double> R;
    vector<double> norms;
    vector<double> exponets;
    vector<double> coefs;
    vector<int> orbital;
    double B;
    double I_A;
    int val_elec;
    int l;
    int m;
    int n;
    double int_pot; // iniozation potential ?? soemthing like that

    /// aaaa need to update based on other map
    // Function()
    //     : atomic_num(0), // Initialize to default value (0 for an integer)
    //       R({}),         // Empty vector
    //       norms({}),
    //       exponets({}),
    //       coefs({}),
    //       orbital({}),
    //       B(0.0),      // Initialize to default value (0.0 for a double)
    //       I_A(0.0),    // Initialize to default value (0.0 for a double)
    //       val_elec(0), // Initialize to default value (0 for an integer)
    //       l(0),        // Initialize to default value (0 for an integer)
    //       m(0),        // Initialize to default value (0 for an integer)
    //       n(0),        // Initialize to default value (0 for an integer)
    //       int_pot(0.0) // Initialize to default value (0.0 for a double)
    // {
    // }

    Function(int atomic_num, vector<double> R, vector<double> norms, vector<double> exponets, vector<double> coefs, vector<int> orbital,
             double B, double I_A, int val_elec, int l, int m, int n, double int_pot)
        : atomic_num(atomic_num), R(R), norms(norms), exponets(exponets), coefs(coefs), orbital(orbital),
          B(B), I_A(I_A), val_elec(val_elec), l(l), m(m), n(n), int_pot(int_pot)
    {
    }
};

/// int is atomic number
// vector<vector<double>> is a vector that first componet will be the coefficants, second will be s, orbitial
// third will be p orbitial if applicable and so on.

// map < int, vector<vector<double, 3>> map;

// vector < pair<int, vector<double>> atoms;
// map<int, Function> atoms;

// basis fucntion centered on U and V

class CNDO2
{
public:
    CNDO2(string file_name);
    // double p_uv(vector<double> C_U, vector<double> C_V, electrons spin);
    // double p_uv_tot(vector<double> C_U, vector<double> C_V);

    double p_uv(arma::Col<double> C_U, arma::Col<double> C_V, electrons spin);
    double p_uv_tot(int U, int V);

    double P_tot(Function &fun);
    // double calc_gamma_AB();
    double gamma_AB(Function &funA, Function &funB);
    double gamma_AA(Function &funA, Function &funB);
    void make_gamma();

    // need to determine valance electrons from the charge ( p and q maybe they are VE)

    void calc_P(arma::Mat<double> &P, arma::Mat<double> &C, electrons spin);
    void calc_f(arma::Mat<double> &f, arma::Mat<double> &P);
    // void calc_C(arma::Mat<double> &C, arma::Mat<double> &f);
    void calc_C_and_e(arma::Mat<double> &C, arma::Col<double> &e, arma::Mat<double> &f);
    void calc_S();
    double inf_norm(arma::Mat<double> &P_new, arma::Mat<double> &P_old);
    void huckle(arma::Mat<double> &H);
    // do not change or code will break
    double total_energy();
    int N;
    pair<double, double> one_step();
    void run();
    void info();

private:
    void read_basis(string filepath, int atomic_num);
    void make_map();
    void add_norms();
    void calc_basis();
    void get_s_basis();
    void read_xyz(string filepath);

    double overlap_integral_1d_analytical(double X_A, double X_B, int l_A, int l_B, double aalpha, double bbeta);
    int binomial_coefficient(int m, int n);
    double X_p(double X_A, double X_B, double aalpha, double bbeta);
    int double_factorial(int num);
    int factorial(int num);

    int alpha;
    int beta;
    vector<pair<int, vector<double>>> atoms;
    vector<Function> basis;
    vector<Function> s_basis;
    // represent the AOs as the rows
    // columns will represent a MO then

    int p;
    int q;
    // double gamma_AB;
    arma::Mat<double> core_H;
    arma::Mat<double> S;
    arma::Mat<double> G;
    arma::Mat<double> F_alpha;
    arma::Mat<double> F_beta;
    arma::Mat<double> P_alpha;
    arma::Mat<double> P_beta;
    arma::Mat<double> C_alpha;
    arma::Mat<double> C_beta;
    arma::Col<double> e_alpha;
    arma::Col<double> e_beta;

    double alpha_change;

    double beta_change;

    int atom_A;
    int atom_B;
    int num_atoms;

    map<int, vector<vector<double>>> basis_map;
    map<int, vector<double>> params_map;
    map<int, int> electron_map;

    // init both still
    // gamma
    double P_AA;
    double P_BB;
};

CNDO2::CNDO2(string file_name)
{

    read_xyz(file_name);
    int a = 0;
    int b = 0;
    int total_val = 0;
    atom_A = 0;
    atom_B = 0;

    for (auto &atom : atoms)
    {
        if (atom.first == 1)
        {
            b++;
            total_val += 1;
        }
        else
        {
            a++;

            // lazy in theroy will only work for elements
            // atomic number 3-10 otherwise need to change code (fine for assignment)
            total_val += (atom.first - 2);
        }
    }

    N = 4 * a + b;

    alpha = total_val / 2;
    beta = alpha;
    if (total_val % 2 != 0)
    {
        alpha++;
    }

    // in order -B, 0.5(I_s + A_s), 0.5(I_p + A_p)
    params_map[1] = {9, 7.176, 0};
    params_map[6] = {21, 14.051, 5.572};
    params_map[7] = {25, 19.316, 7.275};
    params_map[8] = {31, 25.390, 9.111};
    params_map[9] = {39, 32.272, 11.080};

    electron_map[1] = 1;
    electron_map[6] = 4;
    electron_map[7] = 5;
    electron_map[8] = 6;
    electron_map[9] = 7;

    core_H.set_size(N, N);
    S.set_size(N, N);
    F_alpha.set_size(N, N);
    P_alpha.set_size(N, N);
    C_alpha.set_size(N, N);
    F_beta.set_size(N, N);
    P_beta.set_size(N, N);
    C_beta.set_size(N, N);

    // initiate C = zero matrix
    make_map();
    calc_basis();
    add_norms();
    // get_s_basis();

    // G.set_size(s_basis.size(), s_basis.size());
    G.set_size(2, 2);

    calc_S();
    make_gamma();
    huckle(core_H);

    P_alpha.zeros();
    P_beta.zeros();
    alpha_change = 1;
    beta_change = 1;
    // P_AA = P_tot(atom_A);
    // P_BB = P_tot(atom_B);
}

void CNDO2::read_xyz(string filepath)
{

    ifstream myfile(filepath);
    int charge;

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

            // pair<int, vector<double,3>> atom = make_pair(atomic_num, )

            atoms.emplace_back(atomic_num, coords);

            if (atom_A == 0 && atom_B == 0)
            {
                atom_A = atomic_num;
            }
            else if (atom_B == 0)
            {
                atom_B = atomic_num;
            }
        }
    }

    myfile.close();
}

void CNDO2::read_basis(string filepath, int atomic_num)
{
    vector<vector<double>> return_vec;
    vector<double> exponets;

    fstream myfile(filepath);

    if (myfile.is_open())
    {
        double exp;
        double coef;

        if (atomic_num == 1)
        {
            vector<double> coefs;

            for (int i = 0; i < 3; i++)
            {
                myfile >> exp;
                myfile >> coef;
                exponets.push_back(exp);
                coefs.push_back(coef);
            }
            return_vec.push_back(exponets);
            return_vec.push_back(coefs);
        }
        else
        {
            vector<double> coefs_s;
            vector<double> coefs_p;

            for (int i = 0; i < 3; i++)
            {

                myfile >> exp;
                exponets.push_back(exp);
                myfile >> coef;
                coefs_s.push_back(coef);
                myfile >> coef;
                coefs_p.push_back(coef);
            }
            return_vec.push_back(exponets);
            return_vec.push_back(coefs_s);
            return_vec.push_back(coefs_p);
        }
    }

    basis_map[atomic_num] = return_vec;
}

// in order -B, 0.5(I_s + A_s), 0.5(I_p + A_p)
// params_map[1] = {9, 7.176, 0};
// params_map[6] = {21, 14.051, 5.572};
// params_map[7] = {25, 19.316, 7.275};
// params_map[8] = {31, 25.390, 9.111};
// params_map[9] = {39, 32.272, 11.080};

void CNDO2::calc_basis()
{

    for (int i = 0; i < atoms.size(); i++)
    {

        vector<double> norms;
        vector<int> orb = {0, 0, 0};
        ////

        int atom_num = atoms[i].first;
        Function fun = Function(atom_num, atoms[i].second, norms, basis_map[atom_num][0], basis_map[atom_num][1], orb, params_map[atom_num][0], params_map[atom_num][1], electron_map[atom_num], 0, 0, 0, 1);
        basis.push_back(fun);
        s_basis.push_back(fun);
        if (atom_num != 1)
        {
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

                Function fun2 = Function(atom_num, atoms[i].second, norms2, basis_map[atom_num][0], basis_map[atom_num][2], orb2, params_map[atom_num][0], params_map[atom_num][2], electron_map[atom_num], 0, m, 0, 0);

                basis.push_back(fun2);

                count++;
            }
        }
    }
}

double CNDO2::p_uv(arma::Col<double> C_U, arma::Col<double> C_V, electrons spin)
{
    int n;
    arma::Mat<double> temp_C(N, N);

    if (spin == alpha_)
    {
        n = alpha;
    }
    else if (spin == beta_)
    {
        n = beta;
    }

    double sum = 0;
    for (int i = 0; i < n; i++)
    {
        sum += C_U[i] * C_V[i];
    }
    return sum;
}

/// @param C_U should be a row u from the C matrix
/// @param C_V should be a row v from the C matrix
double CNDO2::p_uv_tot(int U, int V)
{
    return p_uv(C_alpha.col(U), C_alpha.col(V), alpha_) + p_uv(C_beta.col(U), C_beta.col(V), beta_);
}

void CNDO2::add_norms()
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

double CNDO2::P_tot(Function &fun)
{
    double sum = 0;
    int u = 0;
    // if (spin == alpha_)
    //{
    for (auto &current_fun : basis)
    {
        // this only works for molecules that are two atom only a
        if (current_fun.R == fun.R)
        {
            sum += p_uv_tot(u, u);
            // sum += P
        }
        u++;
    }
    //}
    // else
    // {
    //     for (auto &fun : basis)
    //     {
    //         if (fun.atomic_num = atomic_num)
    //         {
    //             sum += p_uv_tot(C_beta.col(u), C_beta.col(u));
    //         }
    //         u++;
    //     }
    //}//

    return sum;
}

// will always ve 2x2 for our style problem

// want gamma_XX gamma_XY
//  gamma_XY gamma_YY

// these shoud be Aorbitial of atom_A and Sorbitial

double CNDO2::gamma_AB(Function &funA, Function &funB)
{
    // double RA_RB_squared = pow(funA.R[0] - funB.R[0], 2) + pow(funA.R[1] - funB.R[1], 2) + pow(funA.R[2] - funB.R[2], 2);
    double zero_0;

    double sum = 0;
    for (int k = 0; k < 3; k++)
    {
        for (int kk = 0; kk < 3; kk++)
        {
            for (int l = 0; l < 3; l++)
            {
                for (int ll = 0; ll < 3; ll++)
                {
                    double sigma_a = 1 / (funA.exponets[k] + funA.exponets[kk]);
                    double sigma_b = 1 / (funB.exponets[l] + funB.exponets[ll]);
                    double U_A = pow(M_PI * sigma_a, 1.5);
                    double U_B = pow(M_PI * sigma_b, 1.5);
                    double V_2 = 1 / (sigma_a + sigma_b);
                    double dist = sqrt(pow((funA.R[0] - funB.R[0]), 2) +
                                       pow((funA.R[1] - funB.R[1]), 2) +
                                       pow((funA.R[2] - funB.R[2]), 2));

                    double T = V_2 * pow(dist, 2);
                    zero_0 = U_A * U_B * sqrt(1 / pow(dist, 2)) * erf(sqrt(T));

                    sum += (funA.coefs[k] * funA.norms[k]) * (funA.coefs[kk] * funA.norms[kk]) * (funB.coefs[l] * funB.norms[l]) * (funB.coefs[ll] * funB.norms[ll]) * zero_0;
                }
            }
        }
    }
    return sum * 27.211;
}

void CNDO2::get_s_basis()
{

    for (auto &fun : basis)
    {
        if (fun.orbital == vector<int>{0, 0, 0})
        {
            // Function temp_f = fun;
            s_basis.push_back(fun);
        }
    }
}

double CNDO2::gamma_AA(Function &funA, Function &funB)
{
    // double RA_RB_squared = pow(funA.R[0] - funB.R[0], 2) + pow(funA.R[1] - funB.R[1], 2) + pow(funA.R[2] - funB.R[2], 2);
    double zero_0;

    double sum = 0;
    for (int k = 0; k < 3; k++)
    {
        for (int kk = 0; kk < 3; kk++)
        {
            for (int l = 0; l < 3; l++)
            {
                for (int ll = 0; ll < 3; ll++)
                {
                    double sigma_a = 1 / (funA.exponets[k] + funA.exponets[kk]);
                    double sigma_b = 1 / (funB.exponets[l] + funB.exponets[ll]);
                    double U_A = pow(M_PI * sigma_a, 1.5);
                    double U_B = pow(M_PI * sigma_b, 1.5);
                    double V_2 = 1 / (sigma_a + sigma_b);
                    // double T = V_2 * RA_RB_squared;qu
                    zero_0 = U_A * U_B * sqrt(2 * V_2) * sqrt(2 / M_PI);

                    sum += (funA.coefs[k] * funA.norms[k]) * (funA.coefs[kk] * funA.norms[kk]) * (funB.coefs[l] * funB.norms[l]) * (funB.coefs[ll] * funB.norms[ll]) * zero_0;
                }
            }
        }
    }
    return sum * 27.211;
}

// double CNDO2::gamma_AB(Function &funA, Function &funB)
// {
//     // std::function<double()> zero_0;
//     double zero_0;
//     double sum = 0;
//     double RA_RB_squared = pow(funA.R[0] - funB.R[0], 2) + pow(funA.R[1] - funB.R[1], 2) + pow(funA.R[2] - funB.R[2], 2);

//     for (int k = 0; k < 3; k++)
//     {
//         for (int k_ = 0; k_ < 3; k_++)
//         {
//             for (int l = 0; l < 3; l++)
//             {
//                 for (int l_ = 0; l_ < 3; l_++)
//                 {

//                     if (funA.R != funB.R)
//                     {

//                         // zero_0 = [&funA, &funB, k, k_, l, l_]()

//                         // {

//                         cout
//                             << funA.atomic_num << endl;
//                         for (int i = 0; i < 3; i++)
//                         {
//                             cout << funA.exponets[k] << endl;
//                         }
//                         double sigma_a = 1 / (funA.exponets[k] + funA.exponets[k_]);
//                         double sigma_b = 1 / (funB.exponets[l] + funB.exponets[l_]);
//                         double V_2 = 1 / (sigma_a + sigma_b);

//                         double U_A = pow(M_PI * sigma_a, 3 / 2);
//                         double U_B = pow(M_PI * sigma_b, 3 / 2);
//                         double U = U_A * U_B;

//                         // double RA_RB = pow(pow(funA.R[0] - funB.R[0], 2) + pow(funA.R[1] - funB.R[1], 2) + pow(funA.R[2] - funB.R[2], 2), 0.5);

//                         // double T = V_2 * (pow(funA.R[0] - funB.R[0], 2) + pow(funA.R[1] - funB.R[1], 2) + pow(funA.R[2] - funB.R[2], 2));

//                         // double RA_RB_squared = pow(funA.R[0] - funB.R[0], 2) + pow(funA.R[1] - funB.R[1], 2) + pow(funA.R[2] - funB.R[2], 2);
//                         double T = V_2 * RA_RB_squared;

//                         zero_0 = (U * (1 / sqrt(RA_RB_squared)) * erf(sqrt(T)));
//                         // };
//                     }
//                     else
//                     {
//                         // zero_0 = [&funA, &funB, k, k_, l, l_]()
//                         //{
//                         double sigma_a = 1 / (funA.exponets[k] + funA.exponets[k_]);
//                         double sigma_b = 1 / (funB.exponets[l] + funB.exponets[l_]);
//                         double V_2 = 1 / (sigma_a + sigma_b);

//                         double U_A = pow(M_PI * sigma_a, 3 / 2);
//                         double U_B = pow(M_PI * sigma_b, 3 / 2);
//                         double U = U_A * U_B;

//                         double RA_RB_squared = pow(funA.R[0] - funB.R[0], 2) + pow(funA.R[1] - funB.R[1], 2) + pow(funA.R[2] - funB.R[2], 2);
//                         double T = V_2 * RA_RB_squared;

//                         zero_0 = U * pow(2 * V_2, 0.5) * sqrt(2 / M_PI);
//                         //};
//                     }
//                     sum += (funA.coefs[k] * funA.norms[k]) * (funA.coefs[k_] * funA.norms[k_]) * (funB.coefs[l] * funB.norms[l]) * (funB.coefs[l_] * funB.norms[l_]) * zero_0;
//                 }
//             }
//         }
//     }
//     return sum * 27.211;
// }

// double CNDO2::gamma_AB(){

// }

void CNDO2::make_gamma()
{
    // vector<Function> s_basis;

    // Function & funB;

    // for (Function &fun : basis)
    // {
    //     if (fun.orbital == vector<int>{0, 0, 0})
    //     {
    //         // s_basis.push_back(fun);
    //         s_basis.push_back(fun);
    //     }
    //     // else if (fun.orbital == vector<int>{0, 0, 0} && fun.atomic_num == atom_B)
    //     // {
    //     //     s_basis[1] = fun;
    //     // }
    //     // else
    //     // {
    //     //     continue;
    // }

    // cout << "!!!!!!" << end1;

    // G.at(0, 0) = gamma_AB(basis[0], basis[0]);

    // G.at(0, 1) = gamma_AB(basis[0], basis[1]);
    // G.at(1, 0) = gamma_AB(basis[0], basis[1]);
    // G.at(1, 1) = gamma_AB(basis[1], basis[1]);

    // xx xyx
    // yx yy
    cout << "basis size: " << s_basis.size();
    cout << "N: " << this->N << endl;

    int c1;
    int c2;

    c1 = 0;
    for (auto &fun : basis)
    {
        if (fun.orbital == vector<int>{0, 0, 0})
        {
            // arma::Col<double> temp_col;
            c2 = 0;
            for (auto &fun2 : basis)
            {
                if (c1 == c2)
                {
                    G.at(c1, c2) = gamma_AA(fun, fun);
                }
                else
                {
                    G.at(c1, c2) = gamma_AB(fun, fun2);
                }
                c2++;
            }
            // G.col(c1) = temp_col;
        }
        c1++;
    }

    // for (int i = 0; i < s_basis.size(); i++)
    // {

    //     for (int j = 0; j < s_basis.size(); j++)
    //     {
    //         if (i != j)
    //         {

    //             // cout << "breaking ..." << endl;
    //             // break;
    //             G.at(i, j) = gamma_AB(s_basis[i], s_basis[j]);
    //         }
    //         else
    //         {
    //             G.at(i, j) = gamma_AA(s_basis[i], s_basis[j]);
    //         }
    //     }
    // }

    cout << "!!!" << endl;
    G.print();
}

void CNDO2::calc_P(arma::Mat<double> &P, arma::Mat<double> &C, electrons spin)
{

    int n;
    if (spin == alpha_)
    {
        n = alpha;
    }
    else if (spin == beta_)
    {
        n = beta;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double sum = 0;
            for (int k = 0; k < n; k++)
            {
                sum += C(i, k) * C(j, k);
            }
            P.at(i, j) = sum;
        }
    }
}

void CNDO2::calc_f(arma::Mat<double> &f, arma::Mat<double> &P)
{

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                // f.at(i, j) = 0.5 * (basis[i].B + basis[j].B) * S(i, j) - P(i, j) * gamma_AB(basis[i], basis[j]);
                f.at(i, j) = -0.5 * (basis[i].B + basis[j].B) * S(i, j) - P(i, j) * G(0, 1);
            }
            else
            {

                double sum = 0;
                for (auto &fun : basis)
                {
                    if (fun.R != basis[i].R)
                    {
                        sum += (P_tot(fun) - fun.val_elec) * G(0, 1);
                        // sum += (P_tot(fun) - fun.val_elec) * G() gamma_AB(basis[i], fun);
                    }
                }

                f.at(i, i) = (-1 * (basis[i].I_A) + ((P_tot(basis[i]) - basis[i].val_elec) - (P(i, j) - 0.5)) * G(0, 1)) + sum;
            }

            // should maybe write a vector to keep each atom in individually and then can expand to its basis
            //  i think i might have actually done it this way on assignment 3
            // for (int i = 0; i < num_atoms_in_total; i++)
            //     if (A != B)
            //     {
            //         f.at(i, i) += (P_BB_total - valance_atom_num_B) * gamma_AB;
            //         // here B is any but the current atom we are on
            //         // since this loop is per basis function we need to make sure somehow that this is handled appropiatly

            //         // if element number identifier, as in we can call what atom it is based on the index
            //     }
        }
    }
}

void CNDO2::calc_C_and_e(arma::Mat<double> &C, arma::Col<double> &e, arma::Mat<double> &f)
{

    // arma::Col<double> eigval;
    // arma::Mat<double> eigvec;
    eig_sym(e, C, f);

    // C = eigvec;
}

void CNDO2::make_map()
{
    read_basis("H_STO3G.txt", 1);

    read_basis("C_STO3G.txt", 6);

    read_basis("N_STO3G.txt", 7);

    read_basis("O_STO3G.txt", 8);

    read_basis("F_STO3G.txt", 9);

    // keeps in order but not neccarly, might be good for the first entery tho so we know first map is H.
}

void CNDO2::calc_S()
{
    cout << N << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {

            S.at(i, j) = 0.0;
            Function &u = basis[i];
            Function &v = basis[j];
            cout << "u[i =" << i << "] exp size" << u.exponets.size() << endl;
            cout << "v [j =" << j << "] exp size" << v.exponets.size() << endl;

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

int CNDO2::factorial(int num)
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

int CNDO2::double_factorial(int num)
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

double CNDO2::X_p(double X_A, double X_B, double aalpha, double bbeta)
{
    return (aalpha * X_A + bbeta * X_B) / (aalpha + bbeta);
}

int CNDO2::binomial_coefficient(int m, int n)
{
    return factorial(m) / (factorial(n) * factorial(m - n));
}

double CNDO2::overlap_integral_1d_analytical(double X_A, double X_B, int l_A, int l_B, double aalpha, double bbeta)
{
    double sum = 0;
    double X_P = (aalpha * X_A + bbeta * X_B) / (aalpha + bbeta);
    for (int i = 0; i < l_A + 1; i++)
    {
        for (int j = 0; j < l_B + 1; j++)
        {
            if ((i + j) % 2 == 0)
            {
                sum += binomial_coefficient(l_A, i) * binomial_coefficient(l_B, j) * double_factorial(i + j - 1) * pow(X_P - X_A, l_A - i) * pow(X_P - X_B, l_B - j) / pow(2 * (aalpha + bbeta), (i + j) / 2);
            }
            else
            {
                continue;
            }
        }
    }
    double result = exp(-aalpha * bbeta * pow(X_A - X_B, 2) / (aalpha + bbeta)) * pow(M_PI / (aalpha + bbeta), 0.5) * sum;
    return result;
}

double CNDO2::inf_norm(arma::Mat<double> &P_new, arma::Mat<double> &P_old)
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

void CNDO2::huckle(arma::Mat<double> &H)
{

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i != j)
            {
                // should this be neg??
                H.at(i, j) = 0.5 * ((-1 * basis[i].B) + (-1 * basis[j].B)) * S(i, j);
            }
            else
            {
                int b_num;
                for (auto &atom : atoms)
                {
                    if (atom.first != basis[i].atomic_num)
                    {
                        b_num = atom.first;
                        break;
                    }
                }

                // if (basis[i].atomic_num == 1)
                // {

                //     cout << "atom_A =" << atom_A << endl;
                //     H.at(i, i) = ((-1) * basis[i].I_A) - ((basis[i].val_elec - 0.5) * G(0, 0));
                // }
                // else
                // {
                //     H.at(i, i) = ((-1) * basis[i].I_A) - ((basis[i].val_elec - 0.5) * G(1, 1));
                // }
                H.at(i, i) = ((-1) * basis[i].I_A) - ((basis[i].val_elec - 0.5) * G(0, 0));
                //(electron_map[atoms[1].first] * G(0, 1));

                double sub = 0;
                for (auto &fun : basis)
                {
                    int count1 = 0;
                    int count2 = 0;
                    //&& fun.orbitial == vector<int>{0,0,0}
                    // fun.orbital == vector<int>{0, 0, 0} &&
                    //
                    if (fun.R != basis[i].R)
                    {
                        // sub += fun.val_elec * gamma_AB(basis[i], fun);
                        sub += fun.val_elec * G(0, 1);
                        break;
                    }
                }
                H.at(i, i) -= sub;
            }
        }
    }
}

double CNDO2::total_energy()
{
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            sum += P_alpha(i, j) * (core_H(i, j) + F_alpha(i, j)) + P_beta(i, j) * (core_H(i, j) + F_beta(i, j));
        }
    }
    double dist = pow(pow(atoms[0].second[0] + atoms[1].second[0], 2) + pow(atoms[0].second[1] + atoms[1].second[1], 2) + pow(atoms[0].second[2] + atoms[1].second[2], 2), 0.5);
    sum = 0.5 * sum + ((electron_map[atom_A] * electron_map[atom_B]) / dist);

    return sum;
}

// double CNDO2::inf_norm(arma::Mat<double> &P_new, arma::Mat<double> &P_old)
// {

//     double max_diff = 0;
//     for (int i = 0; i < P_new.n_rows; i++)
//     {
//         for (int j = 0; j < P_new.n_rows; j++)
//         {
//             double diff = abs(P_new(i, j) - P_old(i, j));
//             if (diff > max_diff)
//             {
//                 max_diff = diff;
//             }
//         }
//     }
//     return max_diff;
// }

pair<double, double> CNDO2::one_step()
{

    calc_f(F_alpha, P_alpha);
    calc_f(F_beta, P_beta);

    calc_C_and_e(C_alpha, e_alpha, F_alpha);
    calc_C_and_e(C_beta, e_beta, F_beta);

    arma::Mat<double> new_P_alpha = P_alpha;
    arma::Mat<double> new_P_beta = P_beta;
    calc_P(new_P_alpha, C_alpha, alpha_);
    calc_P(new_P_beta, C_beta, beta_);

    double max_alpha_change = inf_norm(new_P_alpha, P_alpha);
    double max_beta_change = inf_norm(new_P_beta, P_beta);

    this->P_alpha = new_P_alpha;
    this->P_beta = new_P_beta;

    return pair<double, double>(max_alpha_change, max_beta_change);
}

void CNDO2::info()
{

    cout << "gamma:" << endl;
    this->G.print();
    cout << endl;

    cout << "overlap:" << endl;
    S.print();
    cout << endl;

    cout << "core H:" << endl;
    core_H.print();
    cout << endl;
}

void CNDO2::run()
{

    cout << " Starting Guess P_alpha = P_beta = 0" << endl;
    cout << "P_alpha: " << endl;
    P_alpha.print();
    cout << "P_beta: " << endl;
    P_beta.print();

    int itt = 0;
    while (beta_change > 1e-6 || alpha_change > 1e-6)
    {
        cout << endl;
        cout << "iteration: " << itt << endl;
        // cout << "P_alpha: " << endl;
        // this->P_alpha.print();
        // cout << "P_beta: " << endl;
        // this->P_beta.print();

        pair<double, double> max_mag = one_step();

        cout << "F_alpha: " << endl;
        F_alpha.print();
        cout << "F_beta: " << endl;
        F_beta.print();

        cout << "after solving eigen equation: " << itt << endl;

        cout << "C_alpha: " << endl;
        C_alpha.print();
        cout << "C_beta: " << endl;
        C_beta.print();

        cout << "p = " << this->alpha << " and q = " << this->beta << endl;

        cout << "new P_alpha: " << endl;
        P_alpha.print();
        cout << "new P_beta: " << endl;
        P_beta.print();

        alpha_change = max_mag.first;
        beta_change = max_mag.second;
        // pair<double, double> max_mag = one_step();
        itt++;
    }

    cout << "Converged after " << itt << " iterations." << endl;
}

int main()
{

    CNDO2 molecule = CNDO2("H2.txt");
    molecule.info();

    molecule.run();

    double energy = molecule.total_energy();

    cout << endl
         << "Total energy = " << energy << "." << endl;

    CNDO2 molecule2 = CNDO2("HF.txt");
    molecule2.info();

    molecule2.run();

    double energy2 = molecule2.total_energy();

    cout << endl
         << "Total energy = " << energy2 << "." << endl;

    return 0;
}