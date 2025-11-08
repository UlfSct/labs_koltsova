#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

const std::string OUTPUT_DIR = "./results/lab3/";
const std::vector<int> DT_STEPS = { 11, 101, 1001 };
const std::vector<int> DX_STEPS = { 11, 101, 1001 };

const double X_MIN = 0.0;
const double X_MAX = 1.0;
const double T_MIN = 0.0;
const double T_MAX = 1.0;

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

void calculateEquasion(int dt_steps, int dx_steps)
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
        std::vector<double> row;
        row.push_back(7 * n * dt);

        for (int j = 1; j < dx_steps; j++)
        {
            row.push_back((7 * exp(j * dx) * (1 + n * dt) - (f[n][j] - f[n][j - 1]) * ( 1.0f / dx )) * dt + f[n][j]);
        }
        f.push_back(row);
        n++;
    }

    std::ofstream out;
    out.open(OUTPUT_DIR + "n_" + std::to_string(dt_steps) + "_j_" + std::to_string(dx_steps) + ".txt"); 
    printMatrix(f, dt_steps, dx_steps, out);
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
        f.push_back(row);
        f[n + 1][0] = 7 * n * dt;

        for (int j = 1; j < dx_steps - 1; j++)
        {
            f[n + 1][j] = (f[n][j] + dt * (7 * exp(j * dx) * (1 + n * dt)) - dt / (2.0 * dx) * (f[n][j + 1] - f[n + 1][j - 1] - f[n][j])) / (1 + dt / (2.0 * dx));
        }
        int j = dx_steps - 1;
        f[n + 1][j] = (f[n][j] + dt * (7 * exp(j * dx) * (1 + n * dt)) - dt / (2.0 * dx) * (f[n][j] - f[n + 1][j - 1] - f[n][j - 1])) / (1 + dt / (2.0 * dx));
        n++;
    }


    std::ofstream out;
    out.open(OUTPUT_DIR + "imp_n_" + std::to_string(dt_steps) + "_j_" + std::to_string(dx_steps) + ".txt"); 
    printMatrix(f, dt_steps, dx_steps, out);
}


int main()
{
    for (int dt : DT_STEPS)
    {
        for (int dx : DX_STEPS)
        {
            calculateEquasionImplicit(dt, dx);
            if (dx > dt) continue;
            calculateEquasion(dt, dx);
        }
    }

    return 0;
}