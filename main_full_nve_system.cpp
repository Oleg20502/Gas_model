#include <iostream>
#include <fstream>
#include <string>

#include "full_nve_system.h"
#include "saving.h"


template<typename type>
void process(std::string path)
{
    unsigned int N;
    int a;
    type T, dt, tau, x_size, y_size, z_size, eps, sigma, k, rmin, v_mean, v_std;
    bool save;
    std::string par = "Params.txt";
    std::ifstream out(path+par);
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
        k = stod(str1);

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

        getline(out, str2);
        getline(out, str1);
        a = stoi(str1);
    }
    out.close();

    std::cout << N << " particles" << '\n';
    std::cout << T << " model time" << '\n';

    std::ofstream file(path+"System/Parameters_"+std::to_string(a)+".txt");
    if(file.is_open()){
        file << N <<' '<< x_size <<' '<< y_size <<' '<< z_size <<' '<<
            T <<' '<< dt <<' '<< tau <<' '<< eps <<' '<< sigma <<' '<<
            k <<' '<< rmin <<' '<< v_mean <<' '<< v_std <<'\n';
    }
    file.close();

    std::string file1 = "Points/Points_data_"+std::to_string(a)+".txt";
    std::string file2 = "Speed/Speed_data_"+std::to_string(a)+".txt";
    std::string file3 = "System/System_data_"+std::to_string(a)+".txt";

    Full_nve_system s(N, x_size, y_size, z_size, eps, sigma, k);

    //s.set_random_points(rmin);
    s.set_grid_state();
    s.set_random_speed(v_mean, v_std);

    //s.load_points("Data/Points_data.txt");
    //s.load_speed("Data/Speed_data.txt");

    s.init(T, dt, tau);
    s.run(save, path+file1, path+file2, path+file3);

    ChangeStringInFileC(path+par, CountLinesInFile(path+par)-1, std::to_string(a+1));
}


int main()
{
    std::string path = "Full_nve_data/";
    process<double>(path);
    std::cout << "Finished\n";
    return 0;
}

