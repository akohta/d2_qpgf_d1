#include "expint.h"
 
#define MAXIT 100
#define EULER 0.57721566490153286061
#define FPMIN 1.0e-300
#define EPS 1.0e-15


double expint(int n,double x)
{
  int i,ii,nm1;
  double a,b,c,d,del,fact,h,psi,ans;

  nm1=n-1;
  if(n<0 || x<0.0 || (x==0.0 && (n==0 || n==1))){
    printf("bad arguments in expint. n:%d, x:%15.14e\n",n,x);
    return 0.0;
  }
  else {
    if(n==0) return exp(-x)/x;
    else {
      if(x==0.0) return 1.0/nm1;
      else {
        if(x>1.0){
          b=x+n;
          c=1.0/FPMIN;
          d=1.0/b;
          h=d;
          for(i=1;i<=MAXIT;i++){
            a=-(double)(i*(nm1+i));
            b+=2.0;
            d=1.0/(a*d+b);
            c=b+a/c;
            del=c*d;
            h*= del;
            if(fabs(del-1.0)<EPS) return h*exp(-x);
          }
          printf("continued fraction failed in expint\n");
        }
        else {
          ans=(nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
          fact=1.0;
          for(i=1;i<=MAXIT;i++){
            fact*=-x/(double)i;
            if(i!=nm1) del =-fact/(double)(i-nm1);
            else {
              psi=-EULER;
              for(ii=1;ii<=nm1;ii++) psi+=1.0/ii;
              del=fact*(-log(x)+psi);
            }
            ans+=del;
            if(fabs(del)<fabs(ans)*EPS) return ans;
          }
          printf("series failed in expint\n");
        }
      }
    }
  }
  return 0.0;
}

double expint_i(double x)
{
  int k;
  double fact,prev,sum,term;

  if (x <= 0.0) printf("Bad argument in ei. x=%g",x);
  if (x < FPMIN) return log(x)+EULER;
  if (x <= -log(EPS)) {
    sum=0.0;
    fact=1.0;
    for (k=1;k<=MAXIT;k++) {
      fact *= x/k;
      term=fact/k;
      sum += term;
      if (term < EPS*sum) break;
    }
    if (k > MAXIT) printf("Series failed in expint_i()");
    return sum+log(x)+EULER;
  } else {
    sum=0.0;
    term=1.0;
    for (k=1;k<=MAXIT;k++) {
      prev=term;
      term *= k/x;
      if (term < EPS) break;
      if (term < prev) sum += term;
      else {
        sum -= prev;
        break;
      }
    }
    return exp(x)*(1.0+sum)/x;
  }
}
