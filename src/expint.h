#if !defined EXPINT_H
#define EXPINT_H

#include <stdio.h>
#include <math.h>


double expint(int n,double x);
/*
return exponential integral E_n(x)
def E_n(x)=\int_1^infty 1/t^n exp(-xt) dt , x>0 , n=0,1,2,....
*/

double expint_i(double x);
/*
return exponential integral Ei(x),for x>0
def Ei(x)= -P.V.\int_-x^\infty exp(-t)/t dt
*/

#endif
