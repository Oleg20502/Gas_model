#ifndef parameters_h
#define parameters_h

#include <string>

template<typename type>
struct Params{
    unsigned int N = 100;
    type T = 2.0;
    type dt = 0.001;
    type tau = 0.01;

    type x_size = 50.0;
    type y_size = x_size, z_size = x_size;
    //type y_size = 50.0;
    //type z_size = 50.0;

    type eps = 1.0;
    type sigma = 1.0;

    type v_mean = 0.0;
    type v_std = 5.0;

    bool cut = false;

    std::string points_path = "Data/Points_data.txt";
    std::string speed_path = "Data/Speed_data.txt";
    std::string system_path = "Data/System_data.txt";
};

#endif // parameters_h
