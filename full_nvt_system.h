#ifndef full_nvt_system_h
#define full_nvt_system_h

#include "molecule.h"
#include "vector_functions.h"
//#include "space.h"
#include "full_nve_system.h"
#include "LJP.h"

#include <iostream>
#include <vector>
#include <cmath>

using std::sqrt;

class Full_nvt_system: public Full_nve_system {
private:
    int I, J;
    double T, dt, tau;
    double temperature0, xi, Ekin0;

public:
    Full_nvt_system(unsigned int n, double x_size, double y_size, double z_size,
                    double eps, double sigma, double K, double temp, double Xi):
                    Full_nve_system(n, x_size, y_size, z_size, eps, sigma, K)
    {
        temperature0 = temp;
        xi = Xi;
        Ekin0 = 1.5 * N * temperature0;
    }

    inline
    void cor_speed()
    {
        double e = sqrt(1 - dt / xi * (1 - Ekin0/Ekin));
        for(int i=0; i<N; ++i){
            p[i].v_x *= e;
            p[i].v_y *= e;
            p[i].v_z *= e;
        }
    }

    inline
    void count_forces()
    {
        U = 0.0;
        double e = xi*(Ekin0/Ekin - 1);
        for(int i = 0; i<N; ++i){
            F[i][0] = 0;
            F[i][1] = 0;
            F[i][2] = 0;
        }
        for(int i = 0; i<N; ++i){
            for(int j = i+1; j<N; ++j){
                dx = p[i].x - p[j].x;
                if(dx > x_size2) dx -= x_size;
                else if(dx < -x_size2) dx += x_size;

                dy = p[i].y - p[j].y;
                if(dy > y_size2) dy -= y_size;
                else if(dy < -y_size2) dy += y_size;

                dz = p[i].z - p[j].z;
                if(dz > z_size2) dz -= z_size;
                else if(dz < -z_size2) dz += z_size;

                r = sqrt(dx*dx + dy*dy + dz*dz);
                U += LJP(r, eps, sigma);
                Ftmp = LJP_force(r, eps, sigma);
                dFx = Ftmp * dx / r;
                dFy = Ftmp * dy / r;
                dFz = Ftmp * dz / r;
                F[i][0] += dFx;
                F[i][1] += dFy;
                F[i][2] += dFz;

                F[j][0] -= dFx;
                F[j][1] -= dFy;
                F[j][2] -= dFz;
            }
        }
    }

    inline
    void iter(double dt)
    {
        for(int i = 0; i<N; ++i){
            change_speed(i, dt*0.5);
            change_position(i, dt);
        }

        count_forces();

        for(int i = 0; i<N; ++i){
            change_speed(i, dt*0.5);
        }
        cor_speed();
        count_impulse();
        count_kin_energy();
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

        E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
        temperature_e_mean = t = r2 = 0.0;

        for(int i = 0; i<N; ++i){
            pos0[i][0] = p[i].x;
            pos0[i][1] = p[i].y;
            pos0[i][2] = p[i].z;
        }
        pos = pos0;

        count_impulse();
        cor_mean_speed();

        count_kin_energy();
        cor_speed();

        count_impulse();
        count_kin_energy();
        temperature_e = 2.0/3*Ekin/N/k;

        count_forces();
        E = Ekin + U;
    }

    void run(const bool save, std::string file1, std::string file2, std::string file3)
    {
        int counter = 0;
        std::cout << std::scientific;
        std::cout << std::setprecision(7);
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

                //std::cout << "E = " << E_mean << ' ';
                //std::cout << "Impulse = " << Q_mean << ' ';
                std::cout << "Te = " << temperature_e_mean << ' ';
                std::cout << "Ekin = " << Ekin << ' ';
                //std::cout << "rmin = " << find_min(R, N) << ' ';
                //std::cout << "r^2 = " << r2 << ' ';
                std::cout << "Ekin0 =  " << Ekin0 << ' ';
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
