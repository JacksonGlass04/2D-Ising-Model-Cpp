// Import Statements
#include <cmath>
#include <string>
#include <iostream>
#include <time.h>
#include <algorithm>
#include <fstream>

using namespace std;

// Global Constants
const int L = 90;
const double J = 1;
const int N_nn1 = 4;

// Create Lattice Class
class Site
{
public:
    int Sz;

    Site *nn1[N_nn1];
};

// Create Lattice
Site lattice[L][L];

// Set Seed
void init_srand(void)
{
    time_t seconds;
    time(&seconds);
    srand48(seconds);
}

// Set Coordinates
void setCoordinates()
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            lattice[i][j].Sz = 0;
        }
    }
}

// Randomize Spins
void randomizeSpins()
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            double rand = drand48();
            lattice[i][j].Sz = (rand > 0.5) ? 1 : -1;
        }
    }
}

// Mod Function
int mod(int x, int m)
{
    if (x >= 0 && x < m)
        return x;
    else if (x < 0)
        return m - 1 - mod(-1 - x, m);
    else
        return x % m;
}

// Set NN
void set_nn()
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            int h = 0;
            // Up nn
            h = mod((i - 1), L);
            lattice[i][j].nn1[0] = &lattice[h][j];

            // Down nn
            h = mod((i + 1), L);
            lattice[i][j].nn1[1] = &lattice[h][j];

            // Left nn
            h = mod((j + 1), L);
            lattice[i][j].nn1[2] = &lattice[i][h];

            // Right nn
            h = mod((j - 1), L);
            lattice[i][j].nn1[3] = &lattice[i][h];
        }
    }
}

// Total Energy
double Energy()
{
    double sum = 0;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            for (int k = 0; k < N_nn1; k++)
            {
                sum += lattice[i][j].Sz * (lattice[i][j].nn1[k]->Sz);
            }
        }
    }
    return (-J * 0.5 * sum);
}

// Total Magnetization
double Magnetization()
{
    double Magnetization = 0;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            Magnetization += lattice[i][j].Sz;
        }
    }
    return Magnetization;
}

// Boolean Update Spin
bool UpdateSpin(int i, int j, double Temp)
{
    double deltaE = J * (2 * lattice[i][j].Sz) * (lattice[i][j].nn1[0]->Sz + lattice[i][j].nn1[1]->Sz + lattice[i][j].nn1[2]->Sz + lattice[i][j].nn1[3]->Sz);

    double p_rand = drand48();
    if (deltaE == 0)
    {
        if (p_rand <= 0.5)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        double p_flip = exp(-deltaE / Temp);
        if (p_flip > p_rand)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    return false;
}

// Metropolis Algorithm
void MonteCarlo()
{
    // Create Temperature Array
    double T_array[50];

    for (int i = 0; i < 50; i++)
    {
        T_array[i] = 2 + (0.02 * i);
    }

    // nSave/nTherm variables
    int nSave = 500;
    int nTherm = 100000;
    // int nData = 750000;
    int nData = 450000;

    // Open files
    ofstream fs;
    fs.open("E1-L" + to_string(L) + ".txt");
    ofstream eTwo;
    eTwo.open("E2-L" + to_string(L) + ".txt");
    ofstream mOne;
    mOne.open("M1-L" + to_string(L) + ".txt");
    ofstream mTwo;
    mTwo.open("M2-L" + to_string(L) + ".txt");
    ofstream mFour;
    mFour.open("M4-L" + to_string(L) + ".txt");

    // Begin Collecting Data, Loop over T values
    for (int t = 0; t < 50; t++)
    {
        // Output current temp
        cout << "Current temp: ------ " << T_array[t] << "K ------" << endl;
        cout << endl;

        // Randomize spins for each temperature
        setCoordinates();
        randomizeSpins();
        set_nn();

        // Data variables
        double E1 = 0;
        double E2 = 0;
        double M1 = 0;
        double M2 = 0;
        double M4 = 0;

        // Thermal Relaxation
        for (int b = 0; b < nTherm; b++)
        {
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    if (UpdateSpin(i, j, T_array[t]))
                    {
                        lattice[i][j].Sz *= -1;
                    }
                }
            }
        }

        // Sweep over Lattice
        double count = 0;
        for (int b = 0; b < nData; b++)
        {
            // Sweep
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    if (UpdateSpin(i, j, T_array[t]))
                    {
                        lattice[i][j].Sz *= -1;
                    }
                }
            }

            // Every nSave sweeps collect data
            if (b % nSave == 0)
            {
                double energy = Energy();
                E1 = (energy + (E1 * count)) / (count + 1.);
                E2 = (count * E2 + pow(energy, 2)) / (count + 1.);

                double mag = Magnetization();
                M1 = (fabs(mag) + (M1 * count)) / (count + 1.);
                M2 = (count * M2 + pow(mag, 2)) / (count + 1.);
                M4 = (count * M4 + pow(mag, 4)) / (count + 1.);
                count += 1;
            }
        }

        // Write data to each file
        fs << E1 << ',';
        eTwo << E2 << ',';
        mOne << M1 << ',';
        mTwo << M2 << ',';
        mFour << M4 << ',';
    }
    fs.close();
}

// Main Method
int main()
{
    // Set Seed
    init_srand();

    setCoordinates();
    randomizeSpins();

    set_nn();

    MonteCarlo();

    return 0;
}