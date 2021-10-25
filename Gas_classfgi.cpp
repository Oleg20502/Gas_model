#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>

using std::min, std::max, std::sqrt, std::pow, std::abs, std::floor;
using std::stoi, std::stod;

double vec_mod(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

std::string parse_line(std::ifstream out)
{
    std::string str1, str2;
    getline(out, str2);
    getline(out, str1);
    return str1;
}

//Alen Computer Simulation of liquids
double sum(std::vector<double> const &vec)
{
    double sum = 0;
    for(int i = 0; i<vec.size(); ++i)
    {
        sum += vec[i];
    }
    return sum;
}

double find_min(std::vector<std::vector<double>> const &vec, int N)
{
    double min = 1000000000.0;
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            if(vec[i][j] < min && vec[i][j] != 0){
                min = vec[i][j];
            }
        }
    }
    return min;
}

double LJP(double r, double eps=1.0, double sigma=1.0)
{
    double a = pow(sigma/r, 6);
    return 4.0*eps*(a - 1)*a;
}

double LJP_cut(double r, double rmax, double eps=1.0, double sigma=1.0)
{
    if (r < rmax){
        return LJP(r, eps, sigma) - LJP(rmax, eps, sigma);
    }
    else{
        return 0.0;
    }
}

double LJP_force(double r, double eps=1.0, double sigma=1.0)
    {
        double a = pow(sigma/r, 6);
        return 24.0*eps*(2*a - 1)*a/r;
    }


struct Point{
    double m = 1.0;
    double x, y, z;
    double v_x, v_y, v_z;
};


class Space{
private:
    double x_size, y_size, z_size, r, eps = 1.0, sigma = 1.0;
    double E, Ekin, U, Q, Qx, Qy, Qz, P;
    double temperature_e, temperature_m, k = 0.01;
    double E_mean, Ekin_mean, U_mean, Q_mean, Qx_mean, Qy_mean, Qz_mean, P_mean;
    double temperature_e_mean, temperature_m_mean;
    unsigned int N;
    std::vector<Point> p;
    std::vector<Point> px0;
    std::vector<Point> py0;
    std::vector<Point> pz0;
    std::vector<std::vector<double>> R;
    std::vector<std::vector<double>> F;
    std::vector<std::vector<double>> Fx;
    std::vector<std::vector<double>> Fy;
    std::vector<std::vector<double>> Fz;

    double dx, dy, dz;

    std::mt19937_64 rng{time(0)};
    //std::default_random_engine rng{seed};

    double a[100] = {0};


public:
    Space(unsigned int n, double x_size, double y_size, double z_size, double eps, double sigma):
        N{n}, x_size{x_size}, y_size{y_size}, z_size{z_size}, eps{eps}, sigma{sigma},
        p{std::vector<Point> (N)},
        F{std::vector<std::vector<double>> (N, std::vector<double> (N))},
        Fx{std::vector<std::vector<double>> (N, std::vector<double> (N))},
        Fy{std::vector<std::vector<double>> (N, std::vector<double> (N))},
        Fz{std::vector<std::vector<double>> (N, std::vector<double> (N))},
        R{std::vector<std::vector<double>> (N, std::vector<double> (N))}
        {
            for(int i=0; i<N; ++i){
                Fx[i][i] = Fy[i][i] = Fz[i][i] = R[i][i] = F[i][i] = 0.;
            }
        }

    void set_random_points()
    {
        std::uniform_real_distribution <double> disx(0.d , x_size);
        std::uniform_real_distribution <double> disy(0.d , y_size);
        std::uniform_real_distribution <double> disz(0.d , z_size);

        for(int i = 0; i<N; i++)
        {
            p[i].x = disx(rng);
            p[i].y = disy(rng);
            p[i].z = disz(rng);
        }

        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
        }
    }

    void set_crystal_cell()
    {
        if(N==1)
        {
            p[0].x = x_size/2;
            p[0].y = y_size/2;
            p[0].z = z_size/2;
        }

        if(N>1)
        {
            int x = static_cast<int>((log(N)/log(8)))+1;
            int number_of_cubes = pow(8,x);
            float s_x = x_size/pow(2,x);
            float s_y = s_x*y_size/x_size;
            float s_z = s_x*z_size/x_size;
            for(int i = 0; i<pow(2,x); i++)
            {
                for(int j =0; j<pow(2,x); j++)
                {
                    for(int y=0; y<pow(2,x); y++)
                    {
                        p[i*pow(4,x)+j*pow(2,x)+k].x = s_x*(2*i+1);
                        p[i*pow(4,x)+j*pow(2,x)+k].y = s_y*(2*y+1);
                        p[i*pow(4,x)+j*pow(2,x)+k].z = s_z*(2*j+1);

                    }
                }
            }
        }

    }

    void set_grid_points()
    {
        double k = pow(static_cast<double> (0.9*N/x_size/y_size/z_size), 1.0/3);
        int Nx = floor(k*x_size);
        int Ny = floor(k*y_size);
        int Nz = floor(k*z_size);
        std::cout << Nx << '\n';

        for(int i = 0; i<Nx; ++i){
            for(int j=0; j<Ny; ++j){
                for(int m=0; m<Nz; ++m){
                    p[i + j*Nx + m*Nx*Ny].x = i/k+0.1;
                    p[i + j*Nx + m*Nx*Ny].y = j/k+0.1;
                    p[i + j*Nx + m*Nx*Ny].z = m/k+0.1;
                }
            }
        }
        for(int i=N-1; i>Nx*Ny*Nz-1; --i){
            p[i].x = (i+0.5)/k;
            p[i].y = (i+0.5)/k;
            p[i].z = (i+0.5)/k;
        }

        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
        }
    }

    void set_random_speed(double mean, double std)
    {
        std::normal_distribution <double> disv(mean, std);

        for(int i=0; i<N; i++){
            p[i].v_x = disv(rng);
            p[i].v_y = disv(rng);
            p[i].v_z = disv(rng);
        }

        double vc_x, vc_y, vc_z;
        vc_x = vc_y = vc_z = 0.0;
        double M = 0;
        for(int i=0; i<N; i++){
            M += p[i].m;
            vc_x += p[i].v_x * p[i].m;
            vc_y += p[i].v_y * p[i].m;
            vc_z += p[i].v_z * p[i].m;
        }
        vc_x /= M;
        vc_y /= M;
        vc_z /= M;
        for(int i=0; i<N; i++){
            p[i].v_x -= vc_x;
            p[i].v_y -= vc_y;
            p[i].v_z -= vc_z;
        }
    }

    double get_distance(int i, int j)
    {
        dx = abs(p[i].x - p[j].x);
        dy = abs(p[i].y - p[j].y);
        dz = abs(p[i].z - p[j].z);
        return vec_mod(min(dx, x_size-dx), min(dy, y_size-dy), min(dz, z_size-dz));
    }

    void iter(double dt)
    {
        for(int i = 0; i<N; ++i){
            change_speed(i, dt/2);
            change_position(i, dt);
        }

        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
        }

        for(int i = 0; i<N; ++i){
            change_speed(i, dt/2);
        }
    }

    void run(double T, double dt, double tau)
    {
        int I = static_cast<int> (T/dt);
        int J = static_cast<int> (tau/dt);
        int counter = 0;
        E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
        P_mean = temperature_e_mean = temperature_m_mean = 0.0;
        std::cout << std::setprecision(5);
        std::cout << std::fixed;
        for(int i=0; i<I; ++i){
            if (counter == J){
                counter = 0;
                E_mean /= J;
                Ekin_mean /= J;
                U_mean /= J;
                Q_mean /= J;
                Qx_mean /= J;
                Qy_mean /= J;
                Qz_mean /= J;
                temperature_e_mean /= J;
                temperature_m_mean /= J;
                P_mean /= J;
                std::cout << "t = " << i*dt << '\n';
                std::cout << "E = " << E_mean;
                std::cout << "  Impulse = " << Q_mean;
                std::cout << " Te = " << temperature_e_mean;
                std::cout << " Tm = " << temperature_m_mean;
                std::cout << " P = " << P << '\n';
                //std::cout << "  rmin = " << find_min(R, N) << std::endl;
                save_points("Points_data.txt", true);
                E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
                P_mean = temperature_e_mean = temperature_m_mean = 0.0;
            }
            iter(dt);

            count_energy();
            count_impulse();
            count_energy_temperature();
            count_maxwell_temperature();
            count_pressure();

            E_mean += E;
            Ekin_mean += Ekin;
            U_mean += U;
            Q_mean += Q;
            Qx_mean += Qx;
            Qy_mean += Qy;
            Qz_mean += Qz;
            temperature_e_mean += temperature_e;
            temperature_m_mean += temperature_m;
            P_mean += P;
            ++counter;
        }
    }

    void count_forces(int i, int j)
    {
        r = get_distance(i, j);
        R[i][j] = r;
        R[j][i] = r;

        F[i][j] = LJP_force(r, eps, sigma);
//        if (r < 3*eps){
//            F[i][j] = LJP_force(r, eps, sigma);
//        }
//        else{
//            F[i][j] = 0;
//        }
        F[j][i] = F[i][j];

        dx = p[i].x - p[j].x;
        if(dx > x_size/2) dx -= x_size;
        else if(dx < -x_size/2) dx += x_size;

        dy = p[i].y - p[j].y;
        if(dy > y_size/2) dy -= y_size;
        else if(dy < -y_size/2) dy += y_size;

        dz = p[i].z - p[j].z;
        if(dz > z_size/2) dz -= z_size;
        else if(dz < -z_size/2) dz += z_size;

        Fx[i][j] = F[i][j] * dx / r;
        Fy[i][j] = F[i][j] * dy / r;
        Fz[i][j] = F[i][j] * dz / r;

        Fx[j][i] = -Fx[i][j];
        Fy[j][i] = -Fy[i][j];
        Fz[j][i] = -Fz[i][j];
    }

    void change_speed(int i, double dt)
    {
        p[i].v_x += sum(Fx[i])/p[i].m * dt;
        p[i].v_y += sum(Fy[i])/p[i].m * dt;
        p[i].v_z += sum(Fz[i])/p[i].m * dt;
    }

    void change_position(int i, double dt)
    {
        p[i].x += p[i].v_x * dt;
        p[i].y += p[i].v_y * dt;
        p[i].z += p[i].v_z * dt;

        if (p[i].x > x_size) {p[i].x -= x_size;}
        else if (p[i].x < 0) {p[i].x += x_size;}

        if (p[i].y > y_size) {p[i].y -= y_size;}
        else if (p[i].y < 0) {p[i].y += y_size;}

        if (p[i].z > z_size) {p[i].z -= z_size;}
        else if (p[i].z < 0) {p[i].z += z_size;}
    }

    void count_energy()
    {
        Ekin = U = 0.0;
        for(int i = 0; i<N; ++i){
            Ekin += p[i].m * (pow(p[i].v_x, 2) + pow(p[i].v_y, 2) + pow(p[i].v_z, 2)) / 2.0;
            for(int j = i+1; j<N; ++j){
                U += LJP(R[i][j], eps, sigma);
            }
        }
        E = Ekin + U;
    }

    void count_impulse()
    {
        Qx = Qy = Qz = 0.0;
        for(int i=0; i<N; ++i){
            Qx += p[i].v_x * p[i].m;
            Qy += p[i].v_y * p[i].m;
            Qz += p[i].v_z * p[i].m;
        }
        Q = vec_mod(Qx, Qy, Qz);
    }

    void count_energy_temperature()
    {
        temperature_e = 2.0/3*Ekin/N/k;
    }

    void count_maxwell_temperature()
    {
        double v2_msk = 0.0;
        for(int i=0; i<N; ++i){
            v2_msk += pow(p[i].v_x, 2) + pow(p[i].v_y, 2) + pow(p[i].v_z, 2);
        }
        v2_msk /= N;
        temperature_m = v2_msk/3/k;
    }

    void count_pressure()
    {
        P = 0;
        for(int i=0; i<N-1; ++i){
            for(int j=i+1; j<N; ++j){
                P += R[i][j]*F[i][j];
            }
        }
        P /= 3;
        P += N*k*temperature_e_mean/x_size/y_size/z_size;
    }

    void print_points()
    {
        for(int i=0; i<N; ++i){
        std::cout << p[i].x << ' ' << p[i].y << ' ' << p[i].z << ' ' <<
            p[i].v_x << ' ' << p[i].v_y << ' ' << p[i].v_z << '\n';
        }
    }

    void save_points(std::string path, bool save)
    {
        std::ofstream out;
        if (save){
            out.open(path, std::ios::app);
        }
        else{
            out.open(path);
        }
        if (out.is_open()){
            out << N << '\n';
            out << " \n";
            for(int i=0; i<N; ++i){
                out << i << ' ' << p[i].x << ' ' << p[i].y << ' ' << p[i].z << '\n';
            }
        }
        out.close();
    }

    void save_speed(std::string path, bool save)
    {
        std::ofstream out;
        if (save){
            out.open(path, std::ios::app);
        }
        else{
            out.open(path);
        }
        if (out.is_open()){
            //out << N << '\n';
            //out << " \n";
            for(int i=0; i<N; ++i){
                out << i << ' ' << p[i].v_x << ' ' << p[i].v_y << ' ' << p[i].v_z << '\n';
            }
        }
        out.close();
    }

    void load_points(std::string path)
    {
        std::ifstream out(path);
        int a;
        if (out.is_open()){
             out >> N;
             for(int i=0; i<N; ++i){
                out >> a >> p[i].x >> p[i].y >> p[i].z;
             }
        }
    }

    void load_speed(std::string path)
    {
        std::ifstream out(path);
        int a;
        if (out.is_open()){
             out >> N;
             for(int i=0; i<N; ++i){
                out >> a >> p[i].v_x >> p[i].v_y >> p[i].v_z;
             }
        }
    }
};


template<typename type>
void process(std::string path)
{
    unsigned int N;
    type T, dt, tau, x_size, y_size, z_size, eps, sigma, v_mean, v_std;
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
        v_mean = stod(str1);

        getline(out, str2);
        getline(out, str1);
        v_std = stod(str1);
    }

    std::cout << N << " particles" << '\n';

    Space s(N, x_size, y_size, z_size, eps, sigma);
    s.set_random_points();
    //s.set_grid_points();
    s.set_random_speed(v_mean, v_std);

    s.save_points("Points_data.txt", false);

    //s.load_points("Points_data.txt");
    //s.load_speed("Speed_data.txt");

    s.run(T, dt, tau);

    s.save_speed("Speed_data.txt", false);
}


int main()
{
    std::string path = "Params.txt";
    process<double>(path);

    return 0;
}
