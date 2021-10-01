/*
 * d2_qpgf_d1.h
 *
 *  Created on: Oct 17, 2018
 *      Author: ohta
 */

#if !defined D2_QPGF_D1_H
#define D2_QPGF_D1_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>

#include "expint.h"      // Exponential Integral function
#include "Faddeeva.h"    // Faddeeva function
#include "chankel1_01.h" // 1st-kind complex hankel function order 0 and 1  
#include "de_int.h"      // normal DE integrator

// Ewald method param 
#define QG1MAX 100
#define QG2MAX 2500
#define QGELIM 3.0
#define QGCLIM 1.0
#define EXPINTLIMIT 600.0

// struct
typedef struct qp_data{
  double k; // wave number 
  double d; // periodic number ( lattice vector is (d,0) ) 
  double kx; // x component of wave number vector  

}QPDT;


int d2hm_qpgf_d1_ew(double complex *qgf,double complex *dqgf,double *r,double eps,QPDT *td);
/* Ewald method. This is main function.
     qgf : value of quasi-periodic green function.
   *dqgf : value of derivative quasi-periodic green's function.
           dqgf[0]=d(qgf)/dx, dqgf[1]=d(qgf)/dy.
       r : argument (x,y)
           r[0]=x, r[1]=y
     eps : requested relative error.
     *td : pointer of quasi-periodic data ( struct QPDT ).

   return code >0 : normal termination. returned number is total summation number. 
               -1 : function qG1 abnormal termination.
               -2 : function qG2 abnormal termination.
               -3 : out of range ( k <= 0, d <= 0, kx >= k ).
*/

int d2hm_qpgf_d1_ew_cs(double complex *qgf,double complex *dqgf,double *rc,double eps,QPDT *td);
/* Ewald method in cylindrical coordinate system. x=r*cos(theta), y=r*sin(theta)
    *qgf : value of quasi-periodic green function.
   *dqgf : value of derivative quasi-periodic green function.
           dqgf[0]=d(qgf)/dr, dqgf[1]=1/r*d(qgf)/dtheta.
      rc : argument (r,theta).
           rc[0]=r, rc[1]=theta.
     eps : requested relative error.
     *td : pointer of quasi-periodic data ( struct QPDT ).

   return code >0 : normal termination. returned number is total summation number 
               -1 : function qG1 abnormal termination.
               -2 : function qG2 abnormal termination.
               -3 : out of range ( k <= 0, d <= 0, kx >= k ).
*/

int d2hm_qpgf_d1_ew_cs2(double complex *qgf,double complex *dqgf,double *r,double eps,QPDT *td);
/* Ewald method in cylindrical coordinate system. x=r*cos(theta), y=r*sin(theta)
    *qgf : value of quasi-periodic green function.
   *dqgf : value of derivative quasi-periodic green function.
           dqgf[0]=d(qgf)/dr, dqgf[1]=1/r*d(qgf)/dtheta.
       r : argument (x,y). ( orthogonal coordinate )
           r[0]=x, r[1]=y.
     eps : requested relative error.
     *td : pointer of quasi-periodic data ( struct QPDT ).

   return code >0 : normal termination. returned number is total summation number 
               -1 : function qG1 abnormal termination.
               -2 : function qG2 abnormal termination.
               -3 : out of range ( k <= 0, d <= 0, kx >= k ).
*/

int d2hm_qpgf_d1_ir(double complex *qgf,double complex *dqgf,double *r,double eps,QPDT *td,double *err);
/* integral representation. This function will be unstable when k|y| becomes large.  
    *qgf : value of quasi-periodic green function
   *dqgf : value of derivative quasi-periodic green's function
       r : (x,y)
     eps : requested relative error
     *td : pointer of quasi-periodic data
    *err : estimated relative error
           err<0 means abnormal termination of DE integration  
              
   return code >0 : normal termination. returned number is total function calls.
               -1 : function qG abnormal termination.
               -2 : function dqG/dx abnormal termination.
               -3 : function dqG/dy abnormal termination.
               -4 : out of range 
*/

#endif
