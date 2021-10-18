#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>

using std::min, std::max, std::sqrt, std::pow, std::abs, std::floor;

double vec_mod(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

double sim_sum(std::vector<std::vector<double>> const &vec, size_t I)
{
    double sum = 0;
    for(int i = 0; i<I; ++i)
    {
        sum += vec[i][I];
    }
    for(int i = 0; i<vec.size() - I; ++i)
    {
        sum += vec[i][vec.size() - I - 1];
    }
    return sum;
}

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

double force(double r, double eps=1.0, double sigma=1.0)
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
public:
    double x_size, y_size, z_size, r, eps = 1.0, sigma = 1.0;
    unsigned int N;
    std::vector<Point> p;

    std::vector<std::vector<double>> R;
    std::vector<std::vector<double>> F;
    std::vector<std::vector<double>> Fx;
    std::vector<std::vector<double>> Fy;
    std::vector<std::vector<double>> Fz;
    //double k = 1.0;

    double dx, dy, dz;

    std::mt19937_64 rng{time(0)};
    //std::default_random_engine rng{seed};

    double a[100] = {0};


public:
    Space(unsigned int n, double x_size, double y_size, double z_size):
        N{n}, x_size{x_size}, y_size{y_size}, z_size{z_size},
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
                R[i][j] = get_distance(i, j);
                R[j][i] = R[i][j];
            }
        }
    }



    void set_random_speed(double mean, double std)
    {
        std::normal_distribution <double> disv(mean, std);

        for(int i=0; i<N; i++)
        {
            p[i].v_x = disv(rng);
            p[i].v_y = disv(rng);
            p[i].v_z = disv(rng);
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
        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
//            change_speed(i, dt/2);
//            change_position(i, dt);
        }
        for(int i = 0; i<N; ++i){
            change_speed(i, dt/2);
            change_position(i, dt);
        }
        for(int i = 0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                count_forces(i, j);
            }
//            change_speed(i, dt/2);
        }
        for(int i = 0; i<N; ++i){
            change_speed(i, dt/2);
        }
    }

    void run(double T, double dt)
    {
        int I = static_cast<int> (T/dt);
        for(int i=0; i<I; ++i){
            iter(dt);
        }
    }

    void count_forces(int i, int j)
    {
        r = get_distance(i, j);
        R[i][j] = r;
        R[j][i] = r;

        F[i][j] = force(r, eps, sigma);
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

        Fx[j][i] = Fx[i][j];
        Fy[j][i] = Fy[i][j];
        Fz[j][i] = Fz[i][j];
    }

    void change_speed(int i, double dt)
    {
//        p[i].v_x += sim_sum(Fx, i)/p[i].m * dt;
//        p[i].v_y += sim_sum(Fy, i)/p[i].m * dt;
//        p[i].v_z += sim_sum(Fz, i)/p[i].m * dt;

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

    double count_energy()
    {
        double W = 0;
        for(int i = 0; i<N; ++i){
            W = W + p[i].m * (pow(p[i].v_x, 2) + pow(p[i].v_y, 2) + pow(p[i].v_z, 2)) / 2.0;
            for(int j = 0; j<N; ++j){
                if(i!=j){
                    r = R[i][j];
                    W = W + 4*eps*(pow((sigma/r),12) - pow((sigma/r),6)) * p[i].m * p[j].m;
                }
            }
        }
        return W;
    }


    void Maxwell_check()
    {
        std::vector<Point> Vr = p;
        for(int i =0; i<N; i++)
        {
            for(int j = 0; j<N; j++)
            {
                if(vec_mod(Vr[j].v_x, Vr[j].v_y, Vr[j].v_z)<vec_mod(Vr[j-1].v_x, Vr[j-1].v_y, Vr[j-1].v_z))
                {
                    Point a = Vr[j];
                    Vr[j] = Vr[j-1];
                    Vr[j-1] = a;
                }
            }
        }

        double deltav = vec_mod(Vr[N-1].v_x, Vr[N-1].v_y, Vr[N-1].v_z) - vec_mod(Vr[0].v_x, Vr[0].v_y, Vr[0].v_z);



        for(int i=0; i<N; i++)
        {
            a[static_cast< int >(100*(vec_mod(Vr[i].v_x, Vr[i].v_y, Vr[i].v_z)-vec_mod(Vr[0].v_x, Vr[0].v_y, Vr[0].v_z))/deltav)]++;
        }
    }

    void print_points()
    {
        for(int i=0; i<N; ++i){
        std::cout << p[i].x << ' ' << p[i].y << ' ' << p[i].z << ' ' <<
            p[i].v_x << ' ' << p[i].v_y << ' ' << p[i].v_z << '\n';
        }
    }

    void save_points(std::string path)
    {
        std::ofstream out(path);
        if (out.is_open()){
            out << N << '\n';
            out << " \n";
            for(int i=0; i<N; ++i){
                out << i << ' ' << p[i].x << ' ' << p[i].y << ' ' << p[i].z << '\n';
            }
        }
    }

    void save_speed(std::string path)
    {
        std::ofstream out(path);
        if (out.is_open()){
            out << N << '\n';
            out << " \n";
            for(int i=0; i<N; ++i){
                out << i << ' ' << p[i].v_x << ' ' << p[i].v_y << ' ' << p[i].v_z << '\n';
            }
        }
    }

    void load_points(std::string path)
    {
        std::ifstream out(path);

        if (out.is_open()){
             out >> N;
             for(int i=0; i<N; ++i){
                out >> p[i].x >> p[i].y >> p[i].z;
             }
        }
    }

    void load_velocity(std::string path)
    {
        std::ifstream out(path);

        if (out.is_open()){
             out >> N;
             for(int i=0; i<N; ++i){
                out >> p[i].v_x >> p[i].v_y >> p[i].v_z;
             }
        }
    }
};

int main()
{
    int n = 45;
    double T = 1.0;
    double dt = 1/pow(2, 10);
    double tau = 0.01;
    Space s(n, 35.0, 35.0, 35.0);
    s.set_random_points();
    s.set_random_speed(0, 1);

    //std::cout << s.count_energy() << std::endl;

    int I = static_cast<int> (T/dt);
    int J = static_cast<int> (tau/dt);
    for(int i=0; i<I; ++i){
        if (i%J == 0){
            std::cout << "E = " << s.count_energy();
            std::cout << "  rmin = " << find_min(s.R, n) << std::endl;
        }
        s.iter(dt);
    }

    //s.save_points("Points_data.txt");

    return 0;
}
