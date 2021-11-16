#ifndef space_h
#define space_h

#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include "molecule.h"
#include "vector_functions.h"
#include "LJP.h"

using std::min, std::max, std::abs, std::pow, std::floor;

class Space{
protected:
    double x_size, y_size, z_size, r, r2, eps, sigma;
    double x_size2, y_size2, z_size2;
    double dx, dy, dz;
    double E, Ekin, U, Q, Qx, Qy, Qz;
    double temperature_e, k;
    double E_mean, Ekin_mean, U_mean, Q_mean, Qx_mean, Qy_mean, Qz_mean;
    double temperature_e_mean;
    double t;
    double M;
    unsigned int N;
    std::vector<Point<double>> p;
    std::vector<std::vector<double>> pos0;
    std::vector<std::vector<double>> pos;
    std::vector<std::vector<double>> F;
    double Ftmp = 0.0, dFx, dFy, dFz;

    std::mt19937_64 rng{static_cast<long long unsigned int>(time(0))};

public:
    Space(unsigned int n, double x_size, double y_size, double z_size, double eps, double sigma, double K):
        N{n}, x_size{x_size}, y_size{y_size}, z_size{z_size}, eps{eps}, sigma{sigma}, k{K},
        p{std::vector<Point<double>> (N)},
        F{std::vector<std::vector<double>> (N, std::vector<double> (3, 0))},
        pos0{std::vector<std::vector<double>> (N, std::vector<double> (3))}
        {
            x_size2 = x_size/2;
            y_size2 = y_size/2;
            z_size2 = z_size/2;
            M = 0.0;
            for(int i=0; i<N; ++i){
                M+= p[i].m;
            }
        }

    void set_random_points(double rmin)
    {
        std::uniform_real_distribution <double> disx(0.d , x_size);
        std::uniform_real_distribution <double> disy(0.d , y_size);
        std::uniform_real_distribution <double> disz(0.d , z_size);
        int i = 0;
        while(i < N){
            p[i].x = disx(rng);
            p[i].y = disy(rng);
            p[i].z = disz(rng);
            ++i;
            for(int j = 0; j<i-1; ++j){
                if(get_distance(i, j) < rmin){
                    --i;
                    break;
                }
            }
        }
    }

    void set_grid_state()
    {
        int n = floor(pow(N, 1.0/3));
        int nx, ny, nz;
        if(pow(n, 3) == N){
            nx = ny = nz = n;
        }
        else if((n+1)*n*n > N){
            nx = ny = n;
            nz = n + 1;
        }
        else if((n+1)*(n+1)*n > N){
            nx = n;
            ny = nz = n + 1;
        }
        else if(pow(n+1, 3) >= N){
            nx = ny = nz = n + 1;
        }
        dx = x_size/nx;
        dy = y_size/ny;
        dz = z_size/nz;
        int tmp = 0;
        for(int i=0; i<nx; ++i){
            for(int j=0; j<ny; ++j){
                for(int k=0; k<nz; ++k){
                    if(tmp < N){
                        p[tmp].x = dx*(2*i+1)/2;
                        p[tmp].y = dy*(2*j+1)/2;
                        p[tmp].z = dz*(2*k+1)/2;
                    }
                    else{
                        break;
                    }
                    ++tmp;
                }
            }
        }
    }

    void set_random_speed(double mean, double std)
    {
        std::normal_distribution <double> disv(mean, std);
        for(int i=0; i<N; ++i){
            p[i].v_x = disv(rng);
            p[i].v_y = disv(rng);
            p[i].v_z = disv(rng);
        }
    }

    inline
    void cor_mean_speed()
    {
        for(int i=0; i<N; ++i){
            p[i].v_x -= Qx/M;
            p[i].v_y -= Qy/M;
            p[i].v_z -= Qz/M;
        }
    }

    inline
    double get_distance(int i, int j)
    {
        dx = abs(p[i].x - p[j].x);
        dy = abs(p[i].y - p[j].y);
        dz = abs(p[i].z - p[j].z);
        return vec_mod(min(dx, x_size-dx), min(dy, y_size-dy), min(dz, z_size-dz));
    }

    inline
    void change_pot_energy(int i)
    {
        for(int j = i+1; j<N; ++j){
            U += LJP(get_distance(i, j), eps, sigma);
        }
    }

//    inline
//    void count_cut_forces(int i, int j)
//    {
//        r = get_distance(i, j);
//        R[i][j] = r;
//        R[j][i] = r;
//
//        if (r < r_max){
//            F[i][j] = LJP_force(r, eps, sigma);
//            F[j][i] = F[i][j];
//
//            dx = p[i].x - p[j].x;
//            if(dx > x_size/2) dx -= x_size;
//            else if(dx < -x_size/2) dx += x_size;
//
//            dy = p[i].y - p[j].y;
//            if(dy > y_size/2) dy -= y_size;
//            else if(dy < -y_size/2) dy += y_size;
//
//            dz = p[i].z - p[j].z;
//            if(dz > z_size/2) dz -= z_size;
//            else if(dz < -z_size/2) dz += z_size;
//
//            Fx[i][j] = F[i][j] * dx / r;
//            Fy[i][j] = F[i][j] * dy / r;
//            Fz[i][j] = F[i][j] * dz / r;
//
//            Fx[j][i] = -Fx[i][j];
//            Fy[j][i] = -Fy[i][j];
//            Fz[j][i] = -Fz[i][j];
//        }
//        else{
//            F[i][j] = F[j][i] = 0;
//            Fx[i][j] = Fy[i][j] = Fz[i][j] = 0.0;
//            Fx[j][i] = Fy[j][i] = Fz[j][i] = 0.0;
//        }
//    }

    inline
    void change_speed(int i, double dt)
    {
        p[i].v_x += F[i][0]/p[i].m * dt;
        p[i].v_y += F[i][1]/p[i].m * dt;
        p[i].v_z += F[i][2]/p[i].m * dt;
    }

    inline
    void change_position(int i, double dt)
    {
        p[i].x += p[i].v_x * dt;
        p[i].y += p[i].v_y * dt;
        p[i].z += p[i].v_z * dt;

        pos[i][0] += p[i].v_x * dt;
        pos[i][1] += p[i].v_y * dt;
        pos[i][2] += p[i].v_z * dt;

        if (p[i].x > x_size) {p[i].x -= x_size;}
        else if (p[i].x < 0) {p[i].x += x_size;}

        if (p[i].y > y_size) {p[i].y -= y_size;}
        else if (p[i].y < 0) {p[i].y += y_size;}

        if (p[i].z > z_size) {p[i].z -= z_size;}
        else if (p[i].z < 0) {p[i].z += z_size;}
    }

    inline
    void count_kin_energy()
    {
        Ekin = 0.0;
        for(int i=0; i<N; ++i){
            Ekin += p[i].m * (pow(p[i].v_x, 2) + pow(p[i].v_y, 2) + pow(p[i].v_z, 2)) * 0.5;
        }
    }

    inline
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

    void count_r_square_mean()
    {
        r2 = 0.0;
        for(int i = 0; i<N; ++i){
            for(int j = 0; j<3; ++j){
                r2 += pow(pos[i][j] - pos0[i][j], 2);
            }
        }
        r2 /= N;
    }

    void print_points()
    {
        for(int i=0; i<N; ++i){
        std::cout << p[i].x << ' ' << p[i].y << ' ' << p[i].z << ' ' <<
            p[i].v_x << ' ' << p[i].v_y << ' ' << p[i].v_z << '\n';
        }
    }

    void save_points(std::ofstream &out)
    {
        out << N << '\n';
        out << "Lattice=" << x_size << " 0 0 0 " << y_size << " 0 0 0 " << z_size << '\n';
        for(int i=0; i<N; ++i){
            out << i << ' ' << p[i].x << ' ' << p[i].y << ' ' << p[i].z << '\n';
        }
    }

    void save_speed(std::ofstream &out)
    {
        for(int i=0; i<N; ++i){
            out << p[i].v_x << ' ' << p[i].v_y << ' ' << p[i].v_z << '\n';
        }
    }

    void save_system_data(std::ofstream &out)
    {
        out << t << ' ' << E << ' ' << Q << ' ' << temperature_e << ' ' << r2 << '\n';
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
        if (out.is_open()){
             for(int i=0; i<N; ++i){
                out >> p[i].v_x >> p[i].v_y >> p[i].v_z;
             }
        }
    }
};


#endif // space_h
