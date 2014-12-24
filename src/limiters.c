#include "inputs.h"
#include "limiters.h"

REAL minMod(REAL a, REAL b)
{
  return ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0);
}

REAL limiter2(REAL a, REAL b)
{
  #if (LIMITER==MINMOD)
    return minMod(a, b);
  #elif (LIMITER==VAN_LEER)
    if (a*b > 0)
    {
      return (2.*a*b/(a+b));
    } else
    {
      return 0.;
    }
  #elif (LIMITER==MC)
    return minMod(2*minMod(a, b), (a+b)/2.);
  #endif
}

REAL limiter4(REAL a, REAL b, REAL c, REAL d)
{
  return limiter2(limiter2(a, b), limiter2(c, d));
}
