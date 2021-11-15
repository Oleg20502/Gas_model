#ifndef space_h
#define space_h

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <functional>

using std::min, std::max, std::abs, std::pow;

template<typename type>
class Space{
protected:
    type x_size, y_size, z_size, r, r_max, eps, sigma;
    type E, Ekin, U, Q, Qx, Qy, Qz, P, free_path_length;
    type temperature_e, k, d = 1/2, min_R, m = 1.0;
    type E_mean, Ekin_mean, U_mean, Q_mean, Qx_mean, Qy_mean, Qz_mean, P_mean, R_mean;
    type temperature_e_mean;
    type t;
    type r2;
    unsigned int N;
    std::vector<Point<type>> p_old;
    std::vector<Point<type>> p;
    std::vector<std::vector<type>> pos0;
    std::vector<std::vector<type>> pos;
    std::vector<std::vector<type>> R;
    std::vector<std::vector<type>> F;
    std::vector<std::vector<type>> Fx;
    std::vector<std::vector<type>> Fy;
    std::vector<std::vector<type>> Fz;

    //std::function<void(int, int)> count_forces;
    //std::function<void()> count_energy;

    type dx, dy, dz;

    std::mt19937_64 rng{(int)(time(0))};

public:
    Space(unsigned int n, type x_size, type y_size, type z_size, type eps, type sigma, type K):
        N{n}, x_size{x_size}, y_size{y_size}, z_size{z_size}, eps{eps}, sigma{sigma}, k{K},
        p{std::vector<Point<type>> (N)},
        p_old{std::vector<Point<type>> (N)},
        F{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        Fx{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        Fy{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        Fz{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        R{std::vector<std::vector<type>> (N, std::vector<type> (N))},
        pos0{std::vector<std::vector<type>> (N, std::vector<type> (3))}
        {
            for(int i=0; i<N; ++i){
                Fx[i][i] = Fy[i][i] = Fz[i][i] = R[i][i] = F[i][i] = 0.;
            }
        }

    void set_random_points(type rmin)
    {
        std::uniform_real_distribution <type> disx(0.d , x_size);
        std::uniform_real_distribution <type> disy(0.d , y_size);
        std::uniform_real_distribution <type> disz(0.d , z_size);
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


    void set_crystal_cell()
    {
        int n=1;
        int x_d = 0;
        int y_d = 0;
        int z_d = 0;
        while(n<N)
        {
            if(x_d < y_d && x_d < z_d)
            {
                x_d++;
                n = n+(y_d+1)*(z_d+1);
            }
            else if(y_d < x_d && y_d < z_d)
            {
                y_d++;
                n = n+(x_d+1)*(z_d+1);
            }
            else if(z_d < x_d && z_d < y_d)
            {
                z_d++;
                n = n+(y_d+1)*(x_d+1);
            }
            else
            {
                x_d++;
                n = n+(y_d+1)*(z_d+1);
            }
        }
        double dx = x_size/(x_d+1);
        double dy = y_size/(y_d+1);
        double dz = z_size/(z_d+1);

        int tmp = 0;
        for(int i = 0; i<d_x+1; i++){
            for(int j =0; j<d_y+1; j++){
                for(int k=0; k<d_z+1; k++){
                    tmp = i*(d_y+1)*(d_z+1)+j(k+1)+k;
                    if(tmp < N){
                        p[tmp].x = dx*(2*i+1)/2;
                        p[tmp].y = dy*(2*k+1)/2;
                        p[tmp].z = dz*(2*j+1)/2;
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

    inline
    type get_distance(int i, int j)
    {
        dx = abs(p[i].x - p[j].x);
        dy = abs(p[i].y - p[j].y);
        dz = abs(p[i].z - p[j].z);
        return vec_mod(min(dx, x_size-dx), min(dy, y_size-dy), min(dz, z_size-dz));
    }

    inline
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

    inline
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
            Fx[i][j] = Fy[i][j] = Fz[i][j] = 0.0;
            Fx[j][i] = Fy[j][i] = Fz[j][i] = 0.0;
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

    void count_cut_energy()
    {
        Ekin = U = 0.0;
        for(int i = 0; i<N; ++i){
            Ekin += p[i].m * (pow(p[i].v_x, 2) + pow(p[i].v_y, 2) + pow(p[i].v_z, 2)) / 2.0;
            for(int j = i+1; j<N; ++j){
                U += LJP_cut(R[i][j], r_max, eps, sigma);
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

    inline
    void iter(type dt)
    {
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
        }
    }

    void run(type T, type dt, type tau, bool save)
    {
        int I = static_cast<int> (T/dt);
        int J = static_cast<int> (tau/dt);
        int counter = 0;
        r_max = 2.5*sigma;
        E = Ekin = U = Q = P = temperature_e = r2 = 0.0;
        E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
        P_mean = temperature_e_mean = t = 0.0;

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
        }

        std::cout << std::scientific;
        std::cout << std::setprecision(10);
        auto s = std::ios::app;
        if(save){
            s = std::ios::trunc;
        }
        std::ofstream out1("Data/Points_data.txt", s);
        std::ofstream out2("Data/Speed_data.txt", s);
        std::ofstream out3("Data/System_data.txt", s);
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
                P_mean /= J;
                count_r_square_mean();

                //std::cout << "t = " << t << '\n';
                std::cout << "E = " << E_mean << ' ';
                //std::cout << "Impulse = " << Q_mean << ' ';
                std::cout << "Te = " << temperature_e_mean << ' ';
                //std::cout << "P = " << P << ' ';
                //std::cout << "rmin = " << find_min(R, N) << ' ';
                //std::cout << "r^2 = " << r2 << ' ';
                std::cout << t << ' ';
                std::cout << '\n';
                if(save){
                    save_points(out1);
                    save_speed(out2);
                    save_system_data(out3);
                }

                E_mean = Ekin_mean = U_mean = Q_mean = Qx_mean = Qy_mean = Qz_mean = 0.0;
                P_mean = temperature_e_mean = 0.0;
            }

            iter(dt);

            count_energy();
            count_impulse();
            count_energy_temperature();
            //count_pressure();

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
        count_free_path_length_mkt();
        std::cout<<" ";
        //count_free_path_length_r(sum_koll, T);
        //std::cout<<" ";
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
        //out << N << '\n';
        //out << " \n";
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
        //int a;
        if (out.is_open()){
             //out >> N;
             for(int i=0; i<N; ++i){
                //out >> a >> p[i].v_x >> p[i].v_y >> p[i].v_z;
                out >> p[i].v_x >> p[i].v_y >> p[i].v_z;
             }
        }
    }

    type count_min_R()
    {
        min_R = x_size*10;
        for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                if(R[i][j]<min_R)
                {
                    min_R = R[i][j];
                }
            }
        }
        return min_R;
    }

    type count_mean_r()
    {
        type R = 0;
        for(int i=0; i<N; i++)
        {
            for(int j=0; j<N; j++)
            {
                R=R+sqrt((p_old[i].x-p_old[j].x)*(p_old[i].x-p_old[j].x)+(p_old[i].y-p_old[j].y)*(p_old[i].y-p_old[j].y)+(p_old[i].y-p_old[j].y)*(p_old[i].y-p_old[j].y));
            }
        }
        R = R/(N*N);
        return R;
    }

    void count_free_path_length_mkt()
    {
        free_path_length = 1/(sqrt(2)*3.14*(N/(x_size*y_size*z_size))*d*d);
        std::cout<<free_path_length;
    }

    int count_koll_soud()
    {
        int koll = 0;
        //R_mean = count_mean_r;
        for(int i = 0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                if(abs(p[i].v_x - p_old[j].v_x) < abs(p[i].v_x/100) + abs(p_old[j].v_x/100) &&
                abs(p[i].v_y - p_old[j].v_y) < abs(p[i].v_y/100) + abs(p_old[j].v_y/100) &&
                abs(p[i].v_z - p_old[j].v_z) < abs(p[i].v_z/100) + abs(p_old[j].v_z/100) &&
                //sqrt((p_old[i].x-p_old[j].x)*(p_old[i].x-p_old[j].x)+(p_old[i].y-p_old[j].y)*(p_old[i].y-p_old[j].y)+(p_old[i].y-p_old[j].y)*(p_old[i].y-p_old[j].y)) < R_mean/10 &&
                i!=j &&
                R[i][j] == min_R
                )
                {
                    koll++;
                }
            }
        }
        koll = koll/2;
        return koll;
    }

    void count_free_path_length_r(int sum_koll, type T)
    {
        if(sum_koll == 0)
        {
            std::cout<<"Not enough data";
        }
        else
        {
            free_path_length = T*N*(sqrt(2*Ekin/(N*m)))/(sum_koll);
            std::cout<<" "<<free_path_length;
        }
    }
};


#endif // space_h
