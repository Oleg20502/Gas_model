#ifndef space_h
#define space_h

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include <chrono>

using std::min, std::max, std::abs;

template<typename type>
class Space{
private:
    type x_size, y_size, z_size, r, r_max, eps = 1.0, sigma = 1.0;
    type E, Ekin, U, Q, Qx, Qy, Qz, P;
    type temperature_e, k = 0.01;
    type E_mean, Ekin_mean, U_mean, Q_mean, Qx_mean, Qy_mean, Qz_mean, P_mean;
    type temperature_e_mean;
    type t;
    unsigned int N;
    std::vector<Point<type>> p;
    std::vector<Point<type>> px0;
    std::vector<Point<type>> py0;
    std::vector<Point<type>> pz0;
    std::vector<std::vector<type>> R;
    std::vector<std::vector<type>> F;
    std::vector<std::vector<type>> Fx;
    std::vector<std::vector<type>> Fy;
    std::vector<std::vector<type>> Fz;

    type dx, dy, dz;

    std::mt19937_64 rng{time(0)};

public:
    Space(unsigned int n, type x_size, type y_size, type z_size, type eps, type sigma):
        N{n}, x_size{x_size}, y_size{y_size}, z_size{z_size}, eps{eps}, sigma{sigma},
        p{std::vector<Point<type>> (N)},
        F{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        Fx{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        Fy{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        Fz{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        R{std::vector<std::vector<type>> (N, std::vector<type> (N))}
        {
            for(int i=0; i<N; ++i){
                Fx[i][i] = Fy[i][i] = Fz[i][i] = R[i][i] = F[i][i] = 0.;
            }
        }

    void set_random_points()
    {
        std::uniform_real_distribution <type> disx(0.d , x_size);
        std::uniform_real_distribution <type> disy(0.d , y_size);
        std::uniform_real_distribution <type> disz(0.d , z_size);

        for(int i = 0; i<N; i++)
        {
            p[i].x = disx(rng);
            p[i].y = disy(rng);
            p[i].z = disz(rng);
        }
    }

    void set_crystal_cell()
    {
        int n;
        if((log(N)/log(8))-static_cast<int>(log(N)/log(8)) > 0)
        {
            n = static_cast<int>((log(N)/log(8)))+1;
        }
        else
        {
            n = log(N)/log(8);
        }
        int number_of_cubes = pow(8,n);
        float s_x = x_size/pow(2,n);
        float s_y = y_size/pow(2,n);
        float s_z = z_size/pow(2,n);
        int tmp = 0;
        for(int i = 0; i<pow(2,n); i++){
            for(int j =0; j<pow(2,n); j++){
                for(int k=0; k<pow(2,n); k++){
                    tmp = i*pow(4,n)+j*pow(2,n)+k;
                    if(tmp < N){
                        p[tmp].x = s_x*(2*i+1)/2;
                        p[tmp].y = s_y*(2*k+1)/2;
                        p[tmp].z = s_z*(2*j+1)/2;
                    }
                }
            }
        }
    }

    void set_random_speed(type mean, type std)
    {
        std::normal_distribution <type> disv(mean, std);

        for(int i=0; i<N; i++){
            p[i].v_x = disv(rng);
            p[i].v_y = disv(rng);
            p[i].v_z = disv(rng);
        }

        type vc_x, vc_y, vc_z;
        vc_x = vc_y = vc_z = 0.0;
        type M = 0;
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

    type get_distance(int i, int j)
    {
        dx = abs(p[i].x - p[j].x);
        dy = abs(p[i].y - p[j].y);
        dz = abs(p[i].z - p[j].z);
        return vec_mod(min(dx, x_size-dx), min(dy, y_size-dy), min(dz, z_size-dz));
    }

    void iter(type dt)
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

    void run(type T, type dt, type tau)
    {
        int I = static_cast<int> (T/dt);
        int J = static_cast<int> (tau/dt);
        int counter = 0;
        r_max = 2.5*sigma;
        E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
        P_mean = temperature_e_mean = t = 0.0;
        std::cout << std::scientific;
        std::cout << std::setprecision(10);

        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
        }
        save_points("Data/Points_data.txt", false);
        save_speed("Data/Speed_data.txt", false);
        save_system_data("Data/System_data.txt", false);

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
                P_mean /= J;

                //std::cout << "t = " << t << '\n';
                std::cout << "E = " << E_mean;
                std::cout << " Impulse = " << Q_mean;
                std::cout << " Te = " << temperature_e_mean;
                //std::cout << " Tm = " << temperature_m_mean;
                std::cout << " P = " << P << '\n';
                //std::cout << "  rmin = " << find_min(R, N) << std::endl;
                save_points("Data/Points_data.txt", true);
                save_speed("Data/Speed_data.txt", true);
                save_system_data("Data/System_data.txt", true);

                E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
                P_mean = temperature_e_mean = 0.0;
            }
            iter(dt);

            count_energy();
            count_impulse();
            count_energy_temperature();
            count_pressure();

            E_mean += E;
            Ekin_mean += Ekin;
            U_mean += U;
            Q_mean += Q;
            Qx_mean += Qx;
            Qy_mean += Qy;
            Qz_mean += Qz;
            temperature_e_mean += temperature_e;
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

    void count_cut_forces(int i, int j)
    {
        r = get_distance(i, j);
        R[i][j] = r;
        R[j][i] = r;

        if (r < r_max){
            F[i][j] = LJP_force(r, eps, sigma);
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
        else{
            F[i][j] = F[j][i] = 0;
            Fx[i][j] = Fx[i][j] = Fx[i][j] = 0.0;
            Fx[j][i] = Fx[j][i] = Fx[j][i] = 0.0;
        }
    }

    void change_speed(int i, type dt)
    {
        p[i].v_x += sum(Fx[i])/p[i].m * dt;
        p[i].v_y += sum(Fy[i])/p[i].m * dt;
        p[i].v_z += sum(Fz[i])/p[i].m * dt;
    }

    void change_position(int i, type dt)
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
            out << std::setprecision(10);
            out << std::scientific;
            for(int i=0; i<N; ++i){
                out << p[i].v_x << ' ' << p[i].v_y << ' ' << p[i].v_z << '\n';
            }
        }
        out.close();
    }

    void save_system_data(std::string path, bool save)
    {
        std::ofstream out;
        if (save){
            out.open(path, std::ios::app);
        }
        else{
            out.open(path);
        }
        if (out.is_open()){
            out << std::setprecision(10);
            out << std::scientific;
            out << t << ' ' << E << ' ' << Q << ' ' << temperature_e << '\n';
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
        //int a;
        if (out.is_open()){
             //out >> N;
             for(int i=0; i<N; ++i){
                //out >> a >> p[i].v_x >> p[i].v_y >> p[i].v_z;
                out >> p[i].v_x >> p[i].v_y >> p[i].v_z;
             }
        }
    }
};


#endif // space_h
