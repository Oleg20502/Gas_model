#include <cmath>
#include <algorithm>
#include <vector>

using std::min, std::max, std::sqrt, std::pow, std::abs;

double vec_mod(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

double force(double r, double eps=1.0, double sigma=1.0)
    {
        return -24.0*eps*(2*pow(sigma/r, 11) - pow(sigma/r, 5));
    }

struct Point{
    double m = 1.0;
    double x, y, z;
    double v_x, v_y, v_z;
};

class Space{
private:
    double x_size, y_size, z_size;
    unsigned int N;
    std::vector<Point> p;
    //double k = 1.0;

    double dx, dy, dz;
    double r;

public:

    Space(unsigned int n, double x_size, double y_size, double z_size):
           N{n}, x_size{x_size}, y_size{y_size}, z_size{z_size} {}

    void set_points()
    {

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

};

int main()
{
    return 0;
}
