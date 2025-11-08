#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

const std::string OUTPUT_DIR = "./results/lab2/";
const std::vector<int> DT_STEPS = { 11, 101, 1001, 10001 };
const std::vector<int> DX_STEPS = { 2001, 201, 21 };

const double L = 0.2;
const double D = 7 * pow(10, -6);
const double Re = 3;
const double d = 10 * pow(10, -3);

const double X_MIN = 0.0;
const double X_MAX = L;
const double T_MIN = 0.0;
const double T_MAX = 10000.0;

double f_t_0(double x)
{
    return 200 + 50 * x;
}

const double c_s = 10;
const double Pr = 400;
const double c_x_0 = 200;
const double beta = (1.0 + 0.5 * 0.55 * pow(Re, 0.5) * pow(Pr, 1.0 / 3.0)) * D / d;

void printMatrix(std::vector<std::vector<double>> f, int dt_steps, int dx_steps, std::ofstream& out)
{
    out << std::setw(15) << "t\\x";

    for (int j = 0; j < dx_steps; j++)
    {
        out << std::setw(15) << std::fixed << std::setprecision(5) << (X_MAX - X_MIN) / double(dx_steps - 1) * j;
    }
    out << "\n";

    for (int n = dt_steps - 1; n >= 0; n--)
    {
        out << std::setw(15) << std::fixed << std::setprecision(5) << (T_MAX - T_MIN) / double(dt_steps - 1) * n;
        for (int j = 0; j < dx_steps; j++)
        {
            out << std::setw(15) << std::fixed << std::setprecision(5) << f[n][j];
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
        row.push_back(f_t_0(j * dx));
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
        b[0] = c_x_0;

        for (int j = 1; j < dx_steps; j++)
        {
            double a_j = -D * dt / pow(dx, 2);
            double b_j = 1 + 2 * D * dt / pow(dx, 2);
            double c_j = -D * dt / pow(dx, 2);
            double e_j = f[n][j];

            a[j] = -a_j / (b_j + c_j * a[j - 1]);
            b[j] = (e_j - c_j * b[j - 1]) / (b_j + c_j * a[j - 1]);
        }

        f[n + 1][dx_steps - 1] = (b[dx_steps - 1] + dx * beta * c_s / D) / (1 + a[dx_steps - 1] + dx * beta / D);
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