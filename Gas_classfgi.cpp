#include <cmath>
#include <algorithm>
#include <vector>
#include <set>
#include <random>

using std::min, std::max, std::sqrt, std::pow, std::abs;

double vec_mod(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

double force(double r, double eps=1.0, double sigma=1.0)
    {
        return -24.0*eps*(2*pow(sigma/(r*r), 11) - pow(sigma/(r*r), 5));
    }


struct Point{
    double m = 1.0;
    double x, y, z;
    double v_x, v_y, v_z;
};



class Space{
private:
    double x_size, y_size, z_size, r, eps = 1.0, sigma = 1.0;
    unsigned int N;
    std::vector<Point> p;
    //double k = 1.0;

    double dx, dy, dz;

    unsigned seed = 3;

    double a[100] = {0};


public:

    Space(unsigned int n, double x_size, double y_size, double z_size):
           N{n}, x_size{x_size}, y_size{y_size}, z_size{z_size} {

           }

    void set_points()
    {
        std::default_random_engine rng(seed);
        std :: uniform_real_distribution <> disx ( 0.f , x_size ) ;
        std :: uniform_real_distribution <> disy ( 0.f , y_size ) ;
        std :: uniform_real_distribution <> disz ( 0.f , z_size ) ;

        for(int i = 0; i<N; i++)
        {
            p[i].x = disx(rng);
            p[i].y = disy(rng);
            p[i].z = disz(rng);
        }
    }

    void set_v()
    {
       for(int i=0; i<N; i++)
       {
            p[i].v_x = 1;
            p[i].v_y = 1;
            p[i].v_z = 1;
       }
    }

    double get_distance(int i, int j)
    {
        dx = abs(p[i].x - p[j].x);
        dy = abs(p[i].y - p[j].y);
        dz = abs(p[i].z - p[j].z);
        return vec_mod(min(dx, x_size-dx), min(dy, y_size-dy), min(dz, z_size-dz));
    }

    void run(double dt)
    {
        for(int i=0; i<N-1; ++i){
            for(int j = i+1; j<N; ++j){
                r = get_distance(i, j);

                change_speed(i, r, dt/2);
                change_speed(j, r, dt/2);

                change_position(i, dt);
                change_position(j, dt);

                r = get_distance(i, j);

                change_speed(i, r, dt/2);
                change_speed(j, r, dt/2);
            }
        }
    }

    void change_speed(int i, double const &r, double const &dt)
    {
        p[i].v_x += force(r)/p[i].m * p[i].x/r * dt/2;
        p[i].v_y += force(r)/p[i].m * p[i].y/r * dt/2;
        p[i].v_z += force(r)/p[i].m * p[i].z/r * dt/2;
    }

    void change_position(int i, double const &dt)
    {
        p[i].x += p[i].v_x * dt;
        p[i].y += p[i].v_y * dt;
        p[i].z += p[i].v_z * dt;
    }

    double count_energy_of_sysyem()
    {
        double W = 0;
        for(int i = 0; i<N; i++)
        {
            W = W + (p[i].m*vec_mod(p[i].v_x, p[i].v_y, p[i].v_z)*vec_mod(p[i].v_x, p[i].v_y, p[i].v_z ))/2.0;
            for(int j =0; j<N; j++)
            {
                if(i!=j)
                {
                    r = get_distance(i, j);
                    W = W + 4*eps*(pow((sigma/r),12) - pow((sigma/r),6));
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
};

int main()
{
    return 0;
}
