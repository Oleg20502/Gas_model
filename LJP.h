#ifndef LJP_h
#define LJP_h

#include <cmath>

template<typename type>
inline
type LJP(type r, type eps=1.0, type sigma=1.0)
{
    type a = std::pow(sigma/r, 6);
    return 4.0*eps*(a - 1)*a;
}

template<typename type>
inline
type LJP_cut(type r, type rmax, type eps=1.0, type sigma=1.0)
{
    if (r < rmax){
        return LJP(r, eps, sigma) - LJP(rmax, eps, sigma);
    }
    else{
        return (type)0;
    }
}

template<typename type>
inline
type LJP_force(type r, type eps=1.0, type sigma=1.0)
{
    type a = std::pow(sigma/r, 6);
    return 24.0*eps*(2*a - 1)*a/r;
}

#endif /* LJP_h */
