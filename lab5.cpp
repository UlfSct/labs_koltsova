#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

#define M_PI 3.14159265358979323846

const std::string OUTPUT_DIR = "./results/lab5/";
const std::vector<int> DT_STEPS = { 101 };
const std::vector<int> DX_STEPS = { 11 };

const double X_MIN = 0.0;
const double X_MAX = 1.0;
const double T_MIN = 0.0;
const double T_MAX = 1.0;

const double k = 4.4;
const double alpha = 0.5;
const double A = pow(M_PI, alpha) * alpha;

double D(double t)
{
    return 0.02 + t / 5.0;
}

double coeff(double n)
{
    return pow(n + 2, 2) - 4.0 / 3.0;
}

double B(double n, double x, double dt)
{
    return k * (pow(x, 2) * pow(dt, 3.0 / 2.0) * 0.5 * coeff(n) - 2 * D((n + 1) * dt) * pow(n * dt, 2));
}

void printMatrix(std::vector<std::vector<double>> f, int dt_steps, int dx_steps, std::ofstream& out)
{
    out << std::setw(10) << "t\\x";

    for (int j = 0; j < dx_steps; j++)
    {
        out << std::setw(10) << std::fixed << std::setprecision(5) << (X_MAX - X_MIN) / double(dx_steps - 1) * j;
    }
    out << "\n";

    for (int n = dt_steps - 1; n >= 0; n--)
    {
        out << std::setw(10) << std::fixed << std::setprecision(5) << (T_MAX - T_MIN) / double(dt_steps - 1) * n;
        for (int j = 0; j < dx_steps; j++)
        {
            out << std::setw(10) << std::fixed << std::setprecision(5) << f[n][j];
        }
        out << "\n";
    }
    out.close();
}

void calculateEquasionImplicit(int dt_steps, int dx_steps)
{
    double dt = (T_MAX - T_MIN) / double(dt_steps - 1);
    double dx = (X_MAX - X_MIN) / double(dx_steps - 1);
    std::vector<std::vector<double>> f;

    std::vector<double> row;
    for (int j = 0; j < dx_steps; j++)
    {
        row.push_back(0);
    }
    f.push_back(row);

    int n = 0;
    while (n < dt_steps)
    {
        std::vector<double> row(dx_steps, -9999);
        std::vector<double> a(dx_steps, -9999);
        std::vector<double> b(dx_steps, -9999);

        f.push_back(row);
        a[0] = 0;
        b[0] = pow((n + 1) * dt, 2);

        for (int j = 1; j < dx_steps; j++)
        {
            double a_j = -D((n + 1) * dt) * A * pow(dt, alpha) / dx;
            double c_j = -D((n + 1) * dt) * A * pow(dt, alpha) / dx;
            double b_j = 1 + 2 * D((n + 1) * dt) * A * pow(dt, alpha) / dx;
            double e_j = f[n][j] * alpha + B(n, j * dx, dt) * A * pow(dt, alpha);

            a[j] = -a_j / (b_j + c_j * a[j - 1]);
            b[j] = (e_j - c_j * b[j - 1]) / (b_j + c_j * a[j - 1]);
        }

        f[n + 1][dx_steps - 1] = k * pow((n + 1) * dt, 2);
        for (int j = dx_steps - 2; j >= 0; j--)
        {
            f[n + 1][j] = a[j] * f[n + 1][j + 1] + b[j];
        }
        
        n++;
    }

    std::ofstream out;
    out.open(OUTPUT_DIR + "n_" + std::to_string(dt_steps) + "_j_" + std::to_string(dx_steps) + ".txt"); 
    printMatrix(f, dt_steps, dx_steps, out);
}

int main()
{
    for (int dt : DT_STEPS)
    {
        for (int dx : DX_STEPS)
        {
            calculateEquasionImplicit(dt, dx);
        }
    }

    return 0;
}