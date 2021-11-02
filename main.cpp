#include <iostream>
#include <fstream>
#include <string>

#include "molecule.h"
#include "vector_functions.h"
#include "LJP.h"
#include "space.h"


template<typename type>
void process(std::string path)
{
    unsigned int N;
    type T, dt, tau, x_size, y_size, z_size, eps, sigma, rmin, v_mean, v_std;
    bool save;
    std::ifstream out(path);
    std::string str1, str2;
    if (out.is_open()){
        getline(out, str2);
        getline(out, str1);
        N = stoi(str1);

        getline(out, str2);
        getline(out, str1);
        x_size = stod(str1);

        getline(out, str2);
        getline(out, str1);
        y_size = stod(str1);

        getline(out, str2);
        getline(out, str1);
        z_size = stod(str1);

        getline(out, str2);
        getline(out, str1);
        T = stod(str1);

        getline(out, str2);
        getline(out, str1);
        dt = stod(str1);

        getline(out, str2);
        getline(out, str1);
        tau = stod(str1);

        getline(out, str2);
        getline(out, str1);
        eps = stod(str1);

        getline(out, str2);
        getline(out, str1);
        sigma = stod(str1);

        getline(out, str2);
        getline(out, str1);
        rmin = stod(str1);

        getline(out, str2);
        getline(out, str1);
        v_mean = stod(str1);

        getline(out, str2);
        getline(out, str1);
        v_std = stod(str1);

        getline(out, str2);
        getline(out, str1);
        save = stoi(str1);
    }

    std::cout << N << " particles" << '\n';
    std::cout << T << " model time" << '\n';

    Space<double> s(N, x_size, y_size, z_size, eps, sigma);

    //s.set_random_points(rmin);
    s.set_crystal_cell();
    s.set_random_speed(v_mean, v_std);

    //s.load_points("Data/Points_data.txt");
    //s.load_speed("Data/Speed_data.txt");

    s.run(T, dt, tau, save);
}


int main()
{
    std::string path = "Params.txt";
    process<double>(path);
    std::cout << "Finished\n";
    return 0;
}
