#ifndef full_nve_system_h
#define full_nve_system_h

#include "molecule.h"
#include "vector_functions.h"
#include "space.h"
#include "LJP.h"

#include <iostream>
#include <vector>


class Full_nve_system: public Space {
private:
    int I, J;
    double T, dt, tau;

public:
    Full_nve_system(unsigned int n, double x_size, double y_size, double z_size, double eps, double sigma, double K):
                    Space(n, x_size, y_size, z_size, eps, sigma, K) {}

    inline
    void count_forces(int i, int j)
    {
        r = get_distance(i, j);
        R[i][j] = r;
        R[j][i] = r;

        F[i][j] = LJP_force(r, eps, sigma);
        F[j][i] = F[i][j];

        dx = p[i].x - p[j].x;
        if(dx > x_size2) dx -= x_size;
        else if(dx < -x_size2) dx += x_size;

        dy = p[i].y - p[j].y;
        if(dy > y_size2) dy -= y_size;
        else if(dy < -y_size2) dy += y_size;

        dz = p[i].z - p[j].z;
        if(dz > z_size2) dz -= z_size;
        else if(dz < -z_size2) dz += z_size;

        Fx[i][j] = F[i][j] * dx / r;
        Fy[i][j] = F[i][j] * dy / r;
        Fz[i][j] = F[i][j] * dz / r;

        Fx[j][i] = -Fx[i][j];
        Fy[j][i] = -Fy[i][j];
        Fz[j][i] = -Fz[i][j];
    }

    inline
    void change_pot_energy(int i)
    {
        for(int j = i+1; j<N; ++j){
            U += LJP(R[i][j], eps, sigma);
        }
    }

    inline
    void iter(double dt)
    {
        Qx = Qy = Qz = 0.0;
        Ekin = U = 0.0;
        for(int i = 0; i<N; ++i){
            change_speed(i, dt*0.5);
            change_position(i, dt);
        }

        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
        }

        for(int i = 0; i<N; ++i){
            change_speed(i, dt*0.5);
            change_impulse(i);
            change_kin_energy(i);
            change_pot_energy(i);
        }
        Q = vec_mod(Qx, Qy, Qz);
        E = Ekin + U;
        temperature_e = 2.0/3*Ekin/N/k;
    }

    void init(double T0, double dt0, double tau0)
    {
        T = T0;
        dt = dt0;
        tau = tau0;
        I = static_cast<int> (T/dt);
        J = static_cast<int> (tau/dt);

        E = Ekin = U = Q = temperature_e = r2 = 0.0;
        E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
        temperature_e_mean = t = 0.0;

        for(int i = 0; i<N; ++i){
            pos0[i][0] = p[i].x;
            pos0[i][1] = p[i].y;
            pos0[i][2] = p[i].z;
        }
        pos = pos0;

        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
            change_impulse(i);
            change_kin_energy(i);
            change_pot_energy(i);
        }
        Q = vec_mod(Qx, Qy, Qz);
        E = Ekin + U;
        temperature_e = 2.0/3*Ekin/N/k;
    }

    void run(const bool save, std::string file1, std::string file2, std::string file3)
    {
        int counter = 0;
        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        auto s = std::ios::app;
        if(save){
            s = std::ios::trunc;
        }
        std::ofstream out1(file1, s);
        std::ofstream out2(file2, s);
        std::ofstream out3(file3, s);
        if(save){
            out1 << std::setprecision(10);
            out1 << std::scientific;
            out2 << std::setprecision(10);
            out2 << std::scientific;
            out3 << std::setprecision(10);
            out3 << std::scientific;
            save_points(out1);
            save_speed(out2);
            save_system_data(out3);
        }

        for(int i=0; i<I; ++i){
            if (counter == J){
                counter = 0;
                t = i*dt;
                E_mean /= J;
                Ekin_mean /= J;
                U_mean /= J;
                Q_mean /= J;
                Qx_mean /= J;
                Qy_mean /= J;
                Qz_mean /= J;
                temperature_e_mean /= J;
                count_r_square_mean();

                std::cout << "E = " << E_mean << ' ';
                //std::cout << "Impulse = " << Q_mean << ' ';
                std::cout << "Te = " << temperature_e_mean << ' ';
                //std::cout << "rmin = " << find_min(R, N) << ' ';
                //std::cout << "r^2 = " << r2 << ' ';
                std::cout << "time = " << t << ' ';
                std::cout << '\n';
                if(save){
                    save_points(out1);
                    save_speed(out2);
                    save_system_data(out3);
                }

                E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
                temperature_e_mean = 0.0;
            }

            iter(dt);

            E_mean += E;
            Ekin_mean += Ekin;
            U_mean += U;
            Q_mean += Q;
            Qx_mean += Qx;
            Qy_mean += Qy;
            Qz_mean += Qz;
            temperature_e_mean += temperature_e;
            ++counter;
        }
    }
};

#endif // full_nve_system_h
