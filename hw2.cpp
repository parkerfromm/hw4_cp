#include <vector>
#include <iostream>
#include <cmath>
#include <functional>
#include <fstream>
#include <map>
#include <tuple>
using namespace std;

double trapezoid_rule(function<double(double)> f, double a, double b, int n)
{

    double sum = 0.5 * (f(a) + f(b));
    double h = (b - a) / n;
    for (int i = 1; i < n; i++)
    {
        sum += f(a + i * h);
    }

    return sum * h;
}

double overlap_integral_1d_numerical(double X_A, double X_B, double l_A, double l_B, double alpha, double beta)
{
    auto f = [X_A, X_B, l_A, l_B, alpha, beta](double x) -> double
    {
        return pow(x - X_A, l_A) * pow(x - X_B, l_B) * exp(-alpha * pow(x - X_A, 2) - beta * pow(x - X_B, 2));
    };

    double a = 5 * (1 / pow(2 * (alpha * beta) / (alpha + beta), 0.5));

    double result = trapezoid_rule(f, -a, a, 100);
    return result;
}

vector<vector<double>> read_xyz_1(string file_path)
{

    ifstream myfile(file_path);

    vector<vector<double>> result;
    vector<double> gaus;
    double current;
    if (myfile.is_open())
    {
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < 3; i++)
            {
                myfile >> current;
                gaus.push_back(current);
            }
            result.push_back(gaus);
        }
    }
    return result;
}

// double overlap_integral_1d_from_file(void){

// };

int factorial(int num)
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

int double_factorial(int num)
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
double X_p(double X_A, double X_B, double alpha, double beta)
{
    return (alpha * X_A + beta * X_B) / (alpha + beta);
}

int binomial_coefficient(int m, int n)
{
    return factorial(m) / (factorial(n) * factorial(m - n));
}

double overlap_integral_1d_analytical(double X_A, double X_B, int l_A, int l_B, double alpha, double beta)
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

vector<vector<int>> return_orbitals(int angular_momentum)
{
    map<int, vector<vector<int>>> orbitals;
    orbitals[0] = {{0, 0, 0}}; // Wrap 's' in a vector since it's a single configuration
    orbitals[1] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    orbitals[2] = {{2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {1, 1, 0}, {0, 1, 1}, {1, 0, 1}};
    return orbitals[angular_momentum];
}

vector<vector<double>> overlap_matrix(vector<double> X_A, vector<double> X_B, int l_A, int l_B, double alpha, double beta)
{
    vector<vector<double>> matrix;

    vector<vector<int>> X_A_function = return_orbitals(l_A);
    vector<vector<int>> X_B_function = return_orbitals(l_B);
    double result;

    for (int i = 0; i < X_A_function.size(); i++)
    {
        vector<double> row;
        for (int j = 0; j < X_B_function.size(); j++)
        {
            result = 1;
            for (int k = 0; k < 3; k++)
            {
                result *= overlap_integral_1d_analytical(X_A[k], X_B[k], X_A_function[i][k], X_B_function[j][k], alpha, beta);
            }
            row.push_back(result);
        }
        matrix.push_back(row);
    }
    return matrix;
}

vector<tuple<vector<double>, double, int>> read_xyz_2(string file_path)
{

    ifstream myfile(file_path);

    vector<tuple<vector<double>, double, int>> result;
    tuple<vector<double>, int, double> shell;
    double exponet;
    int angular_momentum;
    double current;
    if (myfile.is_open())
    {
        for (int j = 0; j < 2; j++)
        {
            vector<double> coords;
            for (int i = 0; i < 3; i++)
            {
                myfile >> current;
                coords.push_back(current);
            }
            myfile >> exponet;
            myfile >> angular_momentum;
            result.push_back(make_tuple(coords, exponet, angular_momentum));
        }
    }
    return result;
}
void print_2(string file_path)
{
    vector<tuple<vector<double>, double, int>> output = read_xyz_2(file_path);
    vector<double> X_A = get<0>(output[0]);
    double alpha = get<1>(output[0]);
    int l_A = get<2>(output[0]);

    vector<double> X_B = get<0>(output[1]);
    double beta = get<1>(output[1]);
    int l_B = get<2>(output[1]);

    vector<vector<int>> X_A_function = return_orbitals(l_A);
    vector<vector<int>> X_B_function = return_orbitals(l_B);

    vector<vector<double>> matrix = overlap_matrix(X_A, X_B, l_A, l_B, alpha, beta);
    cout << "Shell 1 has " << matrix.size() << " functions." << endl;
    cout << "This shell info: R(" << X_A[0] << ", " << X_A[1] << ", " << X_A[2] << "), with angular momentum: "
         << l_A << ", coefficant: " << alpha << endl;

    cout << endl;
    cout << "Shell 2 has " << matrix[0].size() << " functions." << endl;
    cout << "This shell info: R(" << X_B[0] << ", " << X_B[1] << ", " << X_B[2] << "), with angular momentum: "
         << l_B << ", coefficant: " << beta << endl; // coefficant ??

    cout << endl;

    cout << "Overlap integral between Shell 1 and shell 2" << endl;

    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix[0].size(); j++)
        {
            cout << matrix[i][j] << "     ";
        }
        cout << endl;
    }

    cout << " The componets of angular momentum (l, m, n), for the matrix column, from top to bottom, are listed sequentially as: ";
    for (int i = 0; i < X_A_function.size(); i++)
    {
        cout << " (" << X_A_function[i][0] << ", " << X_A_function[i][1] << ", " << X_A_function[i][2] << ")";
        if (i != X_A_function.size() - 1)
        {
            cout << ",";
        }
    }
    cout << "." << endl;

    cout << " The componets of angular momentum (1, m, n) for the matrix, row from left to right, are listed sequentially as: ";
    for (int i = 0; i < X_B_function.size(); i++)
    {
        cout << " (" << X_B_function[i][0] << ", " << X_B_function[i][1] << ", " << X_B_function[i][2];
        if (i != X_B_function.size() - 1)
        {
            cout << ",";
        }
        cout << "." << endl;
    }
}
int main(void)
{

    // Test case 1: Two s-type functions centered at the origin
    double result1 = overlap_integral_1d_numerical(0, 0, 0, 0, 1, 1);
    cout << "Numerical integral of two s-type functions centered at the origin: " << result1 << endl;

    // Test case 2: An s-type and a p-type function centered at the origin
    double result2 = overlap_integral_1d_numerical(0, 0, 1, 0, 1, 1);
    cout << "Numerical integral of an s-type and a p-type function centered at the origin: " << result2 << endl;

    // Test case 3: Offset by 1
    double offset = 1;
    // Two s-type functions with an offset
    double result3 = overlap_integral_1d_numerical(offset, offset, 0, 0, 1, 1);
    cout << "Numerical integral of two s-type functions with an offset of 1: " << result3 << endl;

    // An s-type and a p-type function with an offset
    double result4 = overlap_integral_1d_numerical(offset, offset, 1, 0, 1, 1);
    cout << "Numerical integral of an s-type and a p-type function with an offset of 1: " << result4 << endl;

    cout << endl;

    double result5 = overlap_integral_1d_analytical(0, 0, 0, 0, 1, 1);
    cout << "Analytical integral of two s-type functions centered at the origin: " << result1 << endl;

    // Test case 2: An s-type and a p-type function centered at the origin
    double result6 = overlap_integral_1d_analytical(0, 0, 0, 1, 1, 1);
    cout << "Analytical integral of an s-type and a p-type function centered at the origin: " << result2 << endl;

    // Test case 3: Offset by 1
    // Two s-type functions with an offset
    double result7 = overlap_integral_1d_analytical(offset, offset, 0, 0, 1, 1);
    cout << "Analytical integral of two s-type functions with an offset of 1: " << result3 << endl;

    // An s-type and a p-type function with an offset
    double result8 = overlap_integral_1d_analytical(offset, offset, 1, 0, 1, 1);
    cout << "Analytical integral of an s-type and a p-type function with an offset of 1: " << result4 << endl;

    // Problem #2

    print_2("2.1.txt");
    cout << endl;
    print_2("2.2.txt");
    cout << endl;
    print_2("2.3.txt");

    // vector<tuple<vector<double>, double, int>> output = read_xyz_2("2.3.txt");
    // vector<double> X_A = get<0>(output[1]);
    // int alpha = get<1>(output[1]);
    // int l_A = get<2>(output[1]);

    // vector<double> X_B = get<0>(output[0]);
    // int beta = get<1>(output[0]);
    // int l_B = get<2>(output[0]);

    // vector<vector<int>> X_A_function = return_orbitals(l_A);
    // vector<vector<int>> X_B_function = return_orbitals(l_B);

    // vector<vector<double>> matrix = overlap_matrix(X_A, X_B, l_A, l_B, alpha, beta);
    // cout << "Shell 1 has " << matrix.size() << " functions." << endl;
    // cout << "This shell infor: R(" << X_A[0] << ", " << X_A[1] << ", " << X_A[2] << "), with angular momentum: "
    //      << alpha << endl; // coefficant ??

    // cout << endl;
    // cout << "Shell 2 has " << matrix[0].size() << " functions." << endl;
    // cout << "This shell infor: R(" << X_B[0] << ", " << X_B[1] << ", " << X_B[2] << "), with angular momentum: "
    //      << beta << endl; // coefficant ??

    // cout << endl;

    // cout << "Overlap integral between Shell 1 and shell 2" << endl;

    // for (int i = 0; i < matrix.size(); i++)
    // {
    //     for (int j = 0; j < matrix[0].size(); j++)
    //     {
    //         cout << matrix[i][j] << "     ";
    //     }
    //     cout << endl;
    // }

    // cout << " The componets of angular momentum (l, m, n), for the matrix column, from top to bottom, are listed sequentially as: ";
    // for (int i = 0; i < X_A_function.size(); i++)
    // {
    //     cout << " (" << X_A_function[i][0] << ", " << X_A_function[i][1] << ", " << X_A_function[i][2] << ")";
    //     if (i != X_A_function.size() - 1)
    //     {
    //         cout << ",";
    //     }
    // }
    // cout << "." << endl;

    // cout << " The componets of angular momentum (1, m, n) for the matrix, row from left to right, are listed sequentially as: ";
    // for (int i = 0; i < X_B_function.size(); i++)
    // {
    //     cout << " (" << X_B_function[i][0] << ", " << X_B_function[i][1] << ", " << X_B_function[i][2];
    //     if (i != X_B_function.size() - 1)
    //     {
    //         cout << ",";
    //     }
    //     cout << "." << endl;
    // }

    return 0;
}
