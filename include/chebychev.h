#ifndef __CHEBYSHEV_H__
#define __CHEBYSHEV_H__
/*
 *	Function calculates Chebyshev Polynomials Tn(x)
 */

namespace Chebyshev
{
  // n = 0
  template <class T> inline T T0(const T& x)
  {
    return static_cast<T>(1.0) ;
  }

  // n = 1
  template <class T> inline T T1(const T& x)
  {
    return x ;
  }

  // n = 2
  template <class T> inline T T2(const T& x)
  {
    return (static_cast<T>(2.0) * x*x) - static_cast<T>(1.0) ;
  }

/*
 *	Tn(x)
 */
  template <class T> inline T Tn(unsigned int n, const T& x)
  {
    if (n == 0)
    {
      return T0<T>(x) ;
    }
    else if (n == 1)
    {
      return T1<T>(x) ;
    }
    else if (n == 2)
    {
      return T2<T>(x) ;
    }

/* We could simply do this:
    return (2.0 * x * Tn(n - 1, x)) - Tn(n - 2, x) ;
   but it could be slow for large n */
 
    T tnm1(T2<T>(x)) ;
    T tnm2(T1<T>(x)) ;
    T tn(tnm1) ;

    for (unsigned int l = 3 ; l <= n ; l++)
    { 
      tn = (static_cast<T>(2.0) * x * tnm1) - tnm2 ;
      tnm2 = tnm1;
      tnm1 = tn;
    }

    return tn ;
  }
}
#endif
