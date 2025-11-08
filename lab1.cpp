#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

const std::string OUTPUT_DIR = "./results/lab1/";
const std::vector<int> DT_STEPS = { 11 };
const std::vector<int> DX_STEPS = { 11 };
const std::vector<int> DY_STEPS = { 31 };

const double X_MIN = 0.0;
const double X_MAX = 1.0;
const double Y_MIN = 0.0;
const double Y_MAX = 3.0;
const double T_MIN = 0.0;
const double T_MAX = 1.0;

const double ksi_0 = 10.0;

const double psi_x_0(double y)
{
    return 2 * y * y;
}

const double psi_x_N(double y)
{
    return 3 + 2 * y * y;
}

const double psi_y_0(double x)
{
    return 3 * x * x;
}

const double psi_y_N(double x)
{
    return 18 + 3 * x * x;
}


const double E = 0.00001;

double norm(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b)
{
    double res = 0.0;
    for (int i = 0; i < a.size(); ++i)
    {
        for (int j = 0; j < a[0].size(); ++j)
        {
            res += pow(a[i][j] - b[i][j], 2);
        }
    }
    return std::sqrt(res);
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

void calculateEquasionImplicit(int dt_steps, int dx_steps, int dy_steps)
{
    double dt = (T_MAX - T_MIN) / double(dt_steps - 1);
    double dx = (X_MAX - X_MIN) / double(dx_steps - 1);
    double dy = (Y_MAX - Y_MIN) / double(dy_steps - 1);

    std::vector<std::vector<std::vector<double>>> ksi;
    std::vector<std::vector<std::vector<double>>> psi;
    psi.push_back(*new std::vector<std::vector<double>>(dx_steps, std::vector<double>(dy_steps, 2)));
    ksi.push_back(*new std::vector<std::vector<double>>(dx_steps, std::vector<double>(dy_steps, ksi_0)));

    // внешний цикл по t
    // дальше в начале считаем одно, потом другое, всё с прогонкой с половинчатым шагом

    int n = 0;
    while (n < dt_steps)
    {
        std::vector<std::vector<double>> ksi_1_2(dx_steps, std::vector<double>(dy_steps, -9999));
        ksi.push_back(*new std::vector<std::vector<double>>(dx_steps, std::vector<double>(dy_steps, -9999)));
        
        for (int j = 0; j < dx_steps; j++)
        {
            for (int k = dy_steps - 1; k >= 0; k--)
            {
                if (j == 0 || k == dy_steps - 1)
                {
                    ksi_1_2[j][k] = ksi_0;
                    continue;
                }
                double v_x = (psi[n][j][k + 1] - psi[n][j][k]) / dy;
                ksi_1_2[j][k] = (ksi[n][j][k] / dt + v_x / dx * ksi_1_2[j - 1][k]) / (1 / dt + v_x / dx);
            }
        }
        for (int j = 0; j < dx_steps; j++)
        {
            for (int k = dy_steps - 1; k >= 0; k--)
            {
                if (j == 0 || k == dy_steps - 1)
                {
                    ksi[n + 1][j][k] = ksi_0;
                    continue;
                }
                double v_y = (psi[n][j][k] - psi[n][j - 1][k]) / dx;
                ksi[n + 1][j][k] = (ksi_1_2[j][k] / dt - v_y / dy * ksi[n + 1][j][k + 1]) / (1 / dt - v_y / dx);
            }
        }  

        std::cout << ksi[n + 1][2][2] << " " << n << " \n";
        
        std::vector<std::vector<double>> psi_0 = psi[n];
        std::vector<std::vector<double>> psi_1_2(dx_steps, std::vector<double>(dy_steps, -9999));
        std::vector<std::vector<double>> psi_1(dx_steps, std::vector<double>(dy_steps, -9999));

        while (norm(psi_0, psi_1) > E)
        {
            if (psi_1[0][0] != -9999) {
                for (int j = 0; j < dx_steps; j++)
                {
                    for (int k = 0; k < dy_steps; k++)
                    {
                        psi_0[j][k] = psi_1[j][k];
                    }
                }
            }

            for (int j = 0; j < dx_steps; j++)
            {
                psi_1_2[j][0] = psi_y_0(j * dx);
            }

            for (int k = 1; k < dy_steps; k++)
            {
                std::vector<double> a(dx_steps, -9999);
                std::vector<double> b(dx_steps, -9999);
                a[0] = 0;
                b[0] = psi_x_0(k * dy);
                for (int j = 1; j < dx_steps; j++)
                {
                    double a_j = dt;
                    double b_j = - (2 * dt + pow(dx, 2));
                    double c_j = dt;
                    double e_j = pow(dx, 2) * (dt * ksi[n][j][k] - psi_0[j][k]);
                    a[j] = -a_j / (b_j + c_j * a[j - 1]);
                    b[j] = (e_j - c_j * b[j - 1]) / (b_j + c_j * a[j - 1]);
                }
                psi_1_2[dx_steps - 1][k] = psi_x_N(k * dy);
                
                for (int j = dx_steps - 2; j >= 0; j--)
                {
                    psi_1_2[j][k] = a[j] * psi_1_2[j + 1][k] + b[j];
                }
            }  
            
            for (int k = 0; k < dy_steps; k++)
            {
                psi_1[0][k] = psi_x_0(k * dy);
            }

            for (int j = 1; j < dx_steps; j++)
            {
                std::vector<double> a(dy_steps, -9999);
                std::vector<double> b(dy_steps, -9999);
                a[0] = 0;
                b[0] = psi_y_0(j * dx);
                for (int k = 1; k < dy_steps; k++)
                {
                    double a_j = dt;
                    double b_j = - (2 * dt + pow(dy, 2));
                    double c_j = dt;
                    double e_j = - pow(dy, 2) * psi_1_2[j][k];
                    a[j] = -a_j / (b_j + c_j * a[j - 1]);
                    b[j] = (e_j - c_j * b[j - 1]) / (b_j + c_j * a[j - 1]);
                }
                psi_1[j][dy_steps - 1] = psi_y_N(j * dx);
            


                for (int k = dy_steps - 2; k >= 0; k--)
                {
                    psi_1[j][k] = a[j] * psi_1_2[j][k + 1] + b[j];
                }
            }
        }

        psi.push_back(*new std::vector<std::vector<double>>(dx_steps, std::vector<double>(dy_steps, -9999)));
        for (int j = 0; j < dx_steps; j++)
        {
            for (int k = 0; k < dy_steps; k++)
            {
                psi[n + 1][j][k] = psi_1[j][k];
            }
        }
        n++;
    }

    std::ofstream out;
    n = 0;
    while (n < dt_steps)
    {
        out.open(OUTPUT_DIR + "ksi_n_" + std::to_string(n) + ".txt"); 
        printMatrix(ksi[n], dx_steps, dy_steps, out);
        out.open(OUTPUT_DIR + "psi_n_" + std::to_string(n) + ".txt"); 
        printMatrix(psi[n], dx_steps, dy_steps, out);
        n++;
    }
    
}

int main()
{
    for (int dt : DT_STEPS)
    {
        for (int dx : DX_STEPS)
        {
            for (int dy : DY_STEPS)
            {
                calculateEquasionImplicit(dt, dx, dy);
            }
        }
    }

    return 0;
}