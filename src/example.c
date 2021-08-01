#include "d2_qpgf_d1.h"

void nsum_qpgf(double complex *qgf,double complex *dqgf,double *r,QPDT *td);

int main()
{
	// init 
	double lambda=1.014; // wave length
	double d=15.1;       // lattice constant ( lattice vector is (d,0) )
	double angle=-0.11;  // parameter for wave number vector
	double eps=1.0e-15;  // requested relative error 

	QPDT td;
	td.k=2.0*M_PI/lambda;  // wave number 
	td.d=d;
	td.kx=td.k*sin(angle); // x component of wave number vector

	double complex ans,qgf,dqgf[2];
	double r[2],err;
	int erc,lm;
	r[0]= 1.11;
	r[1]= 0.12;
  
	// param
	printf("lambda=%g\n",lambda);
	printf("d     =%g\n",d);
	printf("k_x   =%15.14g\n",td.kx);
	printf("r     =(%g,%g)\n",r[0],r[1]);
	printf("eps   =%g\n\n",eps);
  
	// normal sum
	printf("-- normal sum --\n");
	printf("qgf=sum_{l=-inf}^inf G(x+ld,y) exp(il k_x d)\n");
	lm=1000000;
	nsum_qpgf(&qgf,dqgf,r,&td);  
	printf(" qgf   =% 15.14e %+15.14e I, |l|=%d\n",creal(qgf),cimag(qgf),lm);
	printf("dqgf/dx=% 15.14e %+15.14e I\n",creal(dqgf[0]),cimag(dqgf[0]));
	printf("dqgf/dy=% 15.14e %+15.14e I\n",creal(dqgf[1]),cimag(dqgf[1]));	
  
	// Ewald method
	printf(" -- Ewald method --\n");
	erc=d2hm_qpgf_d1_ew(&qgf,dqgf,r,eps,&td);
	printf(" qgf   =% 15.14e %+15.14e I, erc=%d \n",creal(qgf),cimag(qgf),erc);
	printf("dqgf/dx=% 15.14e %+15.14e I\n",creal(dqgf[0]),cimag(dqgf[0]));
	printf("dqgf/dy=% 15.14e %+15.14e I\n",creal(dqgf[1]),cimag(dqgf[1]));

	// integral representation
	printf(" -- integral representation --\n");
	erc=d2hm_qpgf_d1_ir(&qgf,dqgf,r,eps,&td,&err);
	printf(" qgf   =% 15.14e %+15.14e I, erc=%d ,err=%g\n",creal(qgf),cimag(qgf),erc,err);
	printf("dqgf/dx=% 15.14e %+15.14e I\n",creal(dqgf[0]),cimag(dqgf[0]));
	printf("dqgf/dy=% 15.14e %+15.14e I\n",creal(dqgf[1]),cimag(dqgf[1]));

	return 0;
}

void nsum_qpgf(double complex *qgf,double complex *dqgf,double *r,QPDT *td)
{
	double complex tep,epp,epm,h1p,h1m;
	double kxd,krp,krm,y2,xp,xm,srp,srm;
	int l,lm;
  
	lm=1000000;
  
	kxd=td->kx*td->d;
	y2=r[1]*r[1];
	l=0;
	tep=cos(kxd)+I*sin(kxd);
	epp=1.0;
	krp=td->k*sqrt(pow(r[0],2)+y2);
	*qgf=besj0(krp)+I*besy0(krp);
	xp=r[0];
	xm=r[0];
	dqgf[0]=(besj1(krp)+I*besy1(krp))*r[0]/sqrt(r[0]*r[0]+y2);
	dqgf[1]=(besj1(krp)+I*besy1(krp))*r[1]/sqrt(r[0]*r[0]+y2);

	for(l=1;l<lm;l++){
		epp*=tep;
		epm=conj(epp);
		xp+=td->d;
		xm-=td->d;
		srp=sqrt(xp*xp+y2);
		srm=sqrt(xm*xm+y2);
		krp=td->k*srp;
		krm=td->k*srm;
		h1p=besj1(krp)+I*besy1(krp);
		h1m=besj1(krm)+I*besy1(krm);
    
		*qgf+=(besj0(krp)+I*besy0(krp))*epp+(besj0(krm)+I*besy0(krm))*epm;
		dqgf[0]+=h1p*epp*xp/srp+h1m*epm*xp/srm;
		dqgf[1]+=h1p*epp*r[1]/srp+h1m*epm*r[1]/srm;
    
		if(l%(lm/10)==0){
			printf("|l|=%d qgf=% 15.14e %+15.14e I\n",l,creal(I*0.25*(*qgf)),cimag(I*0.25*(*qgf)));
		}
	}
	*qgf*=I*0.25;
	dqgf[0]*=-I*0.25*td->k;
	dqgf[1]*=-I*0.25*td->k;
}


