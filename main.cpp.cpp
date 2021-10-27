#include <iostream>
#include <fstream>

#include "molecule.h"
#include "vector_functions.h"
#include "LJP.h"
#include "space.h"
#include "parameters.h"


int main()
{
    Params<double> par;
    std::cout << par.N << " particles" << '\n';

    Space s(par.N, par.x_size, par.y_size, par.z_size, par.eps, par.sigma);
    //s.set_random_points();
    //s.set_crystal_cell();
    //s.set_random_speed(v_mean, v_std);

    s.load_points("Data/Points_data.txt");
    s.load_speed("Data/Speed_data.txt");

    s.run(par.T, par.dt, par.tau);

    return 0;
}
