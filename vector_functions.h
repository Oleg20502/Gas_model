#ifndef vector_functions_h
#define vector_functions_h

#include <vector>
#include <cmath>

template<typename type>
inline
type vec_mod(type x, type y, type z)
{
    return std::sqrt(x*x + y*y + z*z);
}

template<typename type>
inline
type sum(std::vector<type> const &vec)
{
    type sum = 0.0;
    for(int i = 0; i<vec.size(); ++i)
    {
        sum += vec[i];
    }
    return sum;
}

template<typename type>
inline
type find_min(std::vector<std::vector<type>> const &vec, int N)
{
    type min = 1000000000.0;
    for(int i = 0; i<N; ++i){
        for(int j = 0; j<N; ++j){
            if(vec[i][j] < min && vec[i][j] != 0){
                min = vec[i][j];
            }
        }
    }
    return min;
}

#endif /* vector_functions_h */
