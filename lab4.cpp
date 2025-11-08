#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

const std::string OUTPUT_DIR = "./results/lab4/";
const std::vector<int> DT_STEPS = { 11, 21, 101};
const std::vector<int> DX_STEPS = { 201 };

const double X_MIN = 0.0;
const double X_MAX = 20.0;
const double T_MIN = 0.0;
const double T_MAX = 10.0;

const double L = 11.0;
const double sigma = L / 12.6;
const double m1 = 15 * L / 40;
const double m2 = 25 * L / 40;

double f_t_0(double x)
{
    return exp(-pow((x - m1) / sigma, 2)) + exp(-pow((x - m2) / sigma, 2));
}

double f_x_0()
{
    return f_t_0(0);
}

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
        f.push_back(row);
        f[n + 1][0] = f_x_0();

        for (int j = 1; j < dx_steps - 1; j++)
        {
            f[n + 1][j] = (f[n][j] - dt / (2.0 * dx) * (f[n][j + 1] - f[n + 1][j - 1] - f[n][j])) / (1 + dt / (2.0 * dx));
        }
        int j = dx_steps - 1;
        f[n + 1][j] = (f[n][j] - dt / (2.0 * dx) * (f[n][j] - f[n + 1][j - 1] - f[n][j - 1])) / (1 + dt / (2.0 * dx));
        n++;
    }

    std::ofstream out;
    out.open(OUTPUT_DIR + "z_n_" + std::to_string(dt_steps) + "_j_" + std::to_string(dx_steps) + ".txt"); 
    printMatrix(f, dt_steps, dx_steps, out);
}

void calculateEquasionImplicitLeftCorner(int dt_steps, int dx_steps)
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
        f.push_back(row);
        f[n + 1][0] = f_x_0();

        for (int j = 1; j < dx_steps; j++)
        {
            f[n + 1][j] = (f[n][j] + dt / dx * f[n][j - 1]) / (1 + dt / dx);
        }
        n++;
    }

    std::ofstream out;
    out.open(OUTPUT_DIR + "lc_n_" + std::to_string(dt_steps) + "_j_" + std::to_string(dx_steps) + ".txt"); 
    printMatrix(f, dt_steps, dx_steps, out);
}


int main()
{
    for (int dt : DT_STEPS)
    {
        for (int dx : DX_STEPS)
        {
            calculateEquasionImplicit(dt, dx);
            calculateEquasionImplicitLeftCorner(dt, dx);
        }
    }

    return 0;
}