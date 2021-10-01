#include "d2_qpgf_d1.h"

// struct for de_int
typedef struct hm_qgf_ew_data{
  double c1,c2; 
  int nI[2];    // integer part of cn
  double nF[2]; // fractional part of cn
}QGEA;

typedef struct hm_qgf_ir_data{
  double a,kx,ky,kd;
  double complex eikx[2],cel[2];

  int fc;
}QGDT;


int d2hm_qpgf_d1_ew(double complex *qgf,double complex *dqgf,double *rb,double reps,QPDT *td)
{
  int prm_chk(QPDT *td);
  double veps(double *r,QPDT *td);
  int qgf1_t1(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // normal sum 
  int qgf1_t2(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // DEint sum
  int qgf2_t1(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // normal sum

  double complex tf,tdf[2],cf;
  double ve,kve,r[2];
  int erc1,erc2,ld;

  if(prm_chk(td)){
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -3;
  }
  
  ld=floor((rb[0]+0.5*td->d)/td->d);
  r[0]=rb[0]-(double)ld*td->d;
  r[1]=rb[1];

  ve=veps(r,td);
  kve=td->k*ve;

  // qG1
  if(kve<5.0)  erc1=qgf1_t1(&tf,tdf,r,ve,reps,td); //expint sum
  else         erc1=qgf1_t2(&tf,tdf,r,ve,reps,td); //DEint sum
  if(erc1>=0){
    *qgf=tf;    dqgf[0]=tdf[0];    dqgf[1]=tdf[1];
  }
  else {
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -1;
  }

  // qG2
  erc2=qgf2_t1(&tf,tdf,r,ve,reps,td);
  if(erc2>=0){
    *qgf+=tf;    dqgf[0]+=tdf[0];    dqgf[1]+=tdf[1];
    if(ld==0) return erc1+erc2;
    else {
      cf=cexp(-I*(double)ld*td->kx*td->d);
      *qgf*=cf;
      dqgf[0]*=cf;
      dqgf[1]*=cf;
      return erc1+erc2;
    }
  }
  else {
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -2;
  }
}

int d2hm_qpgf_d1_ew_cs(double complex *qgf,double complex *dqgf,double *rc,double reps,QPDT *td)
{
  int prm_chk(QPDT *td);
  double veps(double *r,QPDT *td);
  int qgf1_t1a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // normal sum 
  int qgf1_t2a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // DEint sum
  int qgf2_t1a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // normal sum

  double complex tf,tdf[2];
  double ve,kve,r[2];
  int erc1,erc2;
  
  if(prm_chk(td)){
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -3;
  }
  
  r[0]=rc[0]*cos(rc[1]);
  r[1]=rc[0]*sin(rc[1]);
  
  ve=veps(r,td);
  kve=td->k*ve;
 
  // qG1
  if(kve<5.0) erc1=qgf1_t1a(&tf,tdf,r,ve,reps,td); //expint sum
  else        erc1=qgf1_t2a(&tf,tdf,r,ve,reps,td); //DEint sum
  if(erc1>=0){
    *qgf=tf;    dqgf[0]=tdf[0];    dqgf[1]=tdf[1];
  }
  else {
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -1;
  }
  
  // qG2
  erc2=qgf2_t1a(&tf,tdf,r,ve,reps,td);
  if(erc2>=0){
    *qgf+=tf;    dqgf[0]+=tdf[0];    dqgf[1]+=tdf[1];
    return erc1+erc2;
  }
  else {
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -2;
  }

}

int d2hm_qpgf_d1_ew_cs2(double complex *qgf,double complex *dqgf,double *r,double reps,QPDT *td)
{
  int prm_chk(QPDT *td);
  double veps(double *r,QPDT *td);
  int qgf1_t1a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // normal sum 
  int qgf1_t2a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // DEint sum
  int qgf2_t1a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td); // normal sum

  double complex tf,tdf[2];
  double ve,kve;
  int erc1,erc2;
  
  if(prm_chk(td)){
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -3;
  }
  
  ve=veps(r,td);
  kve=td->k*ve;
 
  // qG1
  if(kve<5.0) erc1=qgf1_t1a(&tf,tdf,r,ve,reps,td); //expint sum
  else        erc1=qgf1_t2a(&tf,tdf,r,ve,reps,td); //DEint sum
  if(erc1>=0){
    *qgf=tf;    dqgf[0]=tdf[0];    dqgf[1]=tdf[1];
  }
  else {
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -1;
  }
  
  // qG2
  erc2=qgf2_t1a(&tf,tdf,r,ve,reps,td);
  if(erc2>=0){
    *qgf+=tf;    dqgf[0]+=tdf[0];    dqgf[1]+=tdf[1];
    return erc1+erc2;
  }
  else {
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -2;
  }

}

int d2hm_qpgf_d1_ir(double complex *qgf,double complex *dqgf,double *rb,double eps,QPDT *td,double *err)
{
  int prm_chk(QPDT *td);
  int init_qgdt(QGDT *pa,double *r,QPDT *td);
  double complex tfu(double u,void *pa);
  double complex tfx(double u,void *pa);
  double complex tfy(double u,void *pa);

  QGDT pm;
  double complex cf;
  double r0,kr,ier,r[2];
  int ld;

  if(prm_chk(td)){
    *qgf=0.0;    dqgf[0]=0.0;    dqgf[1]=0.0;
    return -4;
  }
  
  ld=floor((rb[0]+0.5*td->d)/td->d);
  r[0]=rb[0]-(double)ld*td->d;
  r[1]=rb[1];

  init_qgdt(&pm,r,td);
  r0=sqrt(r[0]*r[0]+r[1]*r[1]);
  kr=td->k*r0;

  // qG
  *qgf=0.25*I*(besj0(kr)+I*besy0(kr))
     +0.5/M_PI*deintiz(tfu,0.0,&pm,eps,&ier); if(ier<0.0) return -1;
  *err=ier;

  // dqG/dx
  dqgf[0]=-0.25*I*td->k*r[0]*(besj1(kr)+I*besy1(kr))/r0
      +td->k/(2.0*M_PI)*deintiz(tfx,0.0,&pm,eps,&ier); if(ier<0.0) return -2;
  *err+=ier;
  
  // dqG/dy
  dqgf[1]=-0.25*I*td->k*r[1]*(besj1(kr)+I*besy1(kr))/r0
      -td->k/(2.0*M_PI)*deintiz(tfy,0.0,&pm,eps,&ier); if(ier<0.0) return -3;
  *err+=ier;

  *err/=3.0;
  
  if(ld!=0){
    cf=cexp(-I*(double)ld*td->kx*td->d);
    *qgf*=cf;
    dqgf[0]*=cf;
    dqgf[1]*=cf;    
  }
  return pm.fc;
}


///////////////////////////////////////////////////////////////
int prm_chk(QPDT *td)
{
  if(td->k<=0.0) return -1;
  if(td->d<=0.0) return -1;
  if(fabs(td->kx)>=td->k) return -1;
  
  return 0;
}

double veps(double *r,QPDT *td)
{
  double beta_m(QPDT *td);
  
  double ep0,bm;
  double elim=QGELIM;

  ep0=td->d/(2.0*sqrt(M_PI));
  bm=beta_m(td);
  if(td->d>300.0/td->k){
    elim=log(4.0*td->d*bm*QGCLIM);
  }
  if(td->k*td->k*ep0*ep0-r[1]*r[1]/(4.0*ep0*ep0)<=elim) return ep0;
  else {
    return sqrt(elim+sqrt(elim*elim+bm*bm*r[1]*r[1]))/(sqrt(2.0)*bm);
  }
}

double beta_m(QPDT *td) // beta_n maximum
{
  double kd,k2,tmp,an,bn2;
  int nb;

  k2=td->k*td->k;
  kd=2.0*M_PI/td->d;
  tmp=td->kx/kd;
  nb=(int)(tmp<0.0 ? tmp-0.5 : tmp+0.5);

  an=kd*(double)nb-td->kx;
  bn2=k2-an*an;

  return sqrt(bn2);
}

int qgf1_t1(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td)
{
  int s0_qgf1(double *tf1,double *tf0,double kep,double arg,double reps);
  int sl_qgf1(double *af1,double *af0,double kep,double *apm,double reps); 

  double complex tep,epp,dif0,dif1;
  double kep,rp2,rm2,y2,tp1,tp0,af1[2],af0[2],apm[2],i_4e2,xp,xm;
  int l,ef;

  kep=td->k*veps;
  i_4e2=1.0/(4.0*veps*veps);
  y2=r[1]*r[1];
  tep=cos(td->kx*td->d)+I*sin(td->kx*td->d);

  //  l=0;
  epp=1.0;
  rp2=r[0]*r[0]+y2;
  ef=s0_qgf1(&tp1,&tp0,kep,rp2*i_4e2,reps);
  if(ef<0) return -1;
  else{
    *f=epp*tp1;
    df[0]=r[0]*epp*tp0;
    df[1]=epp*tp0;
  }
  if(ef==0) return 0;
  xp=r[0];
  xm=r[0];

  // +-l 
  for(l=1;l<QG1MAX;l++){
    xp+=td->d;
    xm-=td->d;
    epp*=tep;

    rp2=xp*xp+y2;
    rm2=xm*xm+y2;
    apm[0]=rp2*i_4e2;
    apm[1]=rm2*i_4e2;
    ef=sl_qgf1(af1,af0,kep,apm,reps);
    if(ef<0) return -1;
    
    dif0=xp*epp*af0[0]+xm*conj(epp)*af0[1];
    dif1=epp*af0[0]+conj(epp)*af0[1];
    *f +=epp*af1[0]+conj(epp)*af1[1];
    df[0]+=dif0;
    df[1]+=dif1;

    if( (cabs(dif0)==0.0 && cabs(dif1)==0.0) || 
    (cabs(dif0)<cabs(df[0])*reps && cabs(dif1)<cabs(df[1])*reps) ) {
    *f/=4.0*M_PI;
    df[0]*= -1.0/(8.0*M_PI*veps*veps);
    df[1]*=-r[1]/(8.0*M_PI*veps*veps);
    return 2*l+1;
    }
  }
  return -1;
}

int s0_qgf1(double *tf1,double *tf0,double kep,double arg,double reps)
{
  double i_fn,tfn,En1,En0,tkep2,kepn,expa,dif;
  int n;

  expa=exp(-arg);
  i_fn=1.0;
  tkep2=kep*kep;
  kepn=1.0;
  En1=expint(1,arg);
  En0=expint(0,arg);
  if(fabs(En0)<DBL_MIN){
    *tf1=0.0;
    *tf0=0.0;
    return 0;
  }
  *tf1=En1;
  *tf0=En0;
  for(n=1;n<QG1MAX;n++){
    tfn=1.0/(double)n;
    i_fn*=tfn;
    kepn*=tkep2;
    En0=En1;
    En1=tfn*(expa-arg*En0);
    dif  =i_fn*kepn*En0;
    *tf1+=i_fn*kepn*En1;
    *tf0+=dif;
    if(dif==0.0 || fabs(dif)<fabs(*tf0)*reps) return n;
  }
  return -1;
}

int sl_qgf1(double *af1,double *af0,double kep,double *apm,double reps)
{
  double eap,eam,i_fn,tfn,tkep2,kepn,En1p,En1m,En0p,En0m,difp,difm,expintlimit;
  int n;
  
  expintlimit=600.0; // 
  eap=exp(-apm[0]);  eam=exp(-apm[1]);

  i_fn=1.0;
  tkep2=kep*kep;
  kepn=1.0;
  En1p=expint(1,apm[0]);  En1m=expint(1,apm[1]);
  En0p=expint(0,apm[0]);  En0m=expint(0,apm[1]);
  if( (En0p<DBL_MIN && En0m<DBL_MIN) || (apm[0]>expintlimit && apm[1]>expintlimit) ){
    af1[0]=0.0;    af1[1]=0.0;
    af0[0]=0.0;    af0[1]=0.0;
    return 0;
  }
  
  af1[0]=En1p;  af1[1]=En1m;
  af0[0]=En0p;  af0[1]=En0m;
  for(n=1;n<QG1MAX;n++){
    tfn=1.0/(double)n;
    i_fn*=tfn;
    kepn*=tkep2;
    En0p=En1p;
    En1p=tfn*(eap-apm[0]*En0p);
    En0m=En1m;
    En1m=tfn*(eam-apm[1]*En0m);
    difp=i_fn*kepn*En0p;
    difm=i_fn*kepn*En0m;
    af0[0]+=difp;
    af0[1]+=difm;
    af1[0]+=i_fn*kepn*En1p;
    af1[1]+=i_fn*kepn*En1m;
    if( (difp==0.0 && difm==0.0) || (fabs(difp)<=fabs(af0[0])*reps && fabs(difm)<=fabs(af0[1])*reps) ) return n;
  }
  return -1;
}

int qgf1_t2(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td)
{
  double qgf1_inf(double u,void *pa);
  double qgf1_idf(double u,void *pa);

  QGEA tmpd;
  double complex epp,tep,difp,difm,it;
  double rl2,xp,xm,y2,err;
  int l;

  tmpd.c1=td->k*td->k;
  y2=r[1]*r[1];

  // l=0;
  epp=1.0;
  tep=cos(td->kx*td->d)+I*sin(td->kx*td->d);
  rl2=r[0]*r[0]+y2;
  tmpd.c2=0.25*rl2;

  *f=deintd(qgf1_inf,0.0,veps,&tmpd,reps,&err);  if(err<0.0) return -1;
  it=deintd(qgf1_idf,0.0,veps,&tmpd,reps,&err);  if(err<0.0) return -1;
  if( cabs(it)<DBL_MIN ){
    *f*=1.0/(2.0*M_PI);
    df[0]=-r[0]/(4.0*M_PI)*it;
    df[1]=-r[1]/(4.0*M_PI)*it;
    return 0; // under flow
  }
  else {
    df[0]=r[0]*it;
    df[1]=it;
  }
  
  xp=r[0];
  xm=r[0];
  for(l=1;l<QG1MAX;l++){
    xp+=td->d;
    xm-=td->d;
    epp*=tep;
    // +l
    rl2=xp*xp+y2;
    tmpd.c2=0.25*rl2;
    *f +=epp*deintd(qgf1_inf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    difp=epp*deintd(qgf1_idf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    df[0]+=xp*difp;
    // -l
    rl2=xm*xm+y2;
    tmpd.c2=0.25*rl2;
    *f +=conj(epp)*deintd(qgf1_inf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    difm=conj(epp)*deintd(qgf1_idf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    df[0]+=xm*difm;
    df[1]+=difp+difm;
    if( (cabs(difp+difm)==0.0) || cabs(difp+difm)<cabs(df[1])*reps){
    *f*=1.0/(2.0*M_PI);
    df[0]*= -1.0/(4.0*M_PI);
    df[1]*=-r[1]/(4.0*M_PI);
    return 2*l+1;
    }
  }
  return -1;
}

double qgf1_inf(double u,void *pa)
{
  QGEA *td=(QGEA *)pa;
  double te,i_u=1.0/u;
  te=exp(td->c1*u*u-td->c2*i_u*i_u);
  return i_u*te;
}

double qgf1_idf(double u,void *pa)
{
  QGEA *td=(QGEA *)pa;
  double te,i_u=1.0/u;
  te=exp(td->c1*u*u-td->c2*i_u*i_u);
  return i_u*i_u*i_u*te;
}

int qgf2_t1(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td)
{
  int sr_qgf2(double complex *tnf,double complex *tyf,double bn,double ay,double veps,double reps);
  int si_qgf2(double complex *tnf,double complex *tyf,double cn,double ay,double veps,double reps);
  int calc_cn(QPDT *td,QGEA *ta);
  double calc_bn2(int n,QGEA *ta);

  QGEA ta;
  double complex dif1,dif2,eiax,tnf=0.0,tyf=0.0,tc;
  double k2,kd,anp,anm,bn,cn,tmp,kl,kh;
  int nb,n;
  
  kl=0.8*td->k;
  kh=1.2*td->k;

  k2=td->k*td->k;
  kd=2.0*M_PI/td->d;
  tmp=td->kx/kd;
  nb=(int)(tmp<0.0 ? tmp-0.5 : tmp+0.5);

  // init QGEA
  if(calc_cn(td,&ta)!=0) return -2;
  ta.c1=4.0*M_PI/(td->d*td->d);
  ta.c2=td->k*td->d;

  // n=nb;
  anp=kd*(double)nb-td->kx;
  anm=anp;
  if(fabs(anp)<kl || fabs(anp)>kh)   tmp=k2-anp*anp;
  else tmp=calc_bn2(nb,&ta);
  eiax=cos(anp*r[0])+I*sin(anp*r[0]);
  if(tmp>0.0){
    bn=sqrt(tmp);
    sr_qgf2(&tnf,&tyf,bn,fabs(r[1]),veps,reps);
    tc=eiax/bn*tnf;
    *f=tc;
    df[0]=anp*tc;
    df[1]=eiax*tyf;
  }
  else if(tmp<0.0){
    cn=sqrt(-tmp);
    si_qgf2(&tnf,&tyf,cn,fabs(r[1]),veps,reps);
    tc=eiax/(I*cn)*tnf;
    *f=tc;
    df[0]=anp*tc;
    df[1]=eiax*tyf;
  }
  else return -2; // bn==0.0

  // nb+-n                                                                                                                                                                        
  for(n=1;n<QG2MAX;n++){
    anp+=kd;
    anm-=kd;
    // +
    if(fabs(anp)<kl || fabs(anp)>kh)   tmp=k2-anp*anp;
    else tmp=calc_bn2(nb+n,&ta);
    eiax=cos(anp*r[0])+I*sin(anp*r[0]);
    if(tmp>0.0){
      bn=sqrt(tmp);
      sr_qgf2(&tnf,&tyf,bn,fabs(r[1]),veps,reps);
      tc=eiax/bn*tnf;
      dif1=anp*tc;
      dif2=eiax*tyf;
      *f+=tc;
    }
    else if(tmp<0.0){
      cn=sqrt(-tmp);
      si_qgf2(&tnf,&tyf,cn,fabs(r[1]),veps,reps);
      tc=eiax/(I*cn)*tnf;
      dif1=anp*tc;
      dif2=eiax*tyf;
      *f+=tc;
    }
    else return -2;
    
    // -
    if(fabs(anm)<kl || fabs(anm)>kh)   tmp=k2-anm*anm;
    else tmp=calc_bn2(nb-n,&ta);
    eiax=cos(anm*r[0])+I*sin(anm*r[0]);
    if(tmp>0.0){
      bn=sqrt(tmp);
      sr_qgf2(&tnf,&tyf,bn,fabs(r[1]),veps,reps);
      tc=eiax/bn*tnf;
      dif1+=anm*tc;
      dif2+=eiax*tyf;
      *f+=tc;
    }
    else if(tmp<0.0){
      cn=sqrt(-tmp);
      si_qgf2(&tnf,&tyf,cn,fabs(r[1]),veps,reps);
      tc=eiax/(I*cn)*tnf;
      dif1+=anm*tc;
      dif2+=eiax*tyf;
      *f+=tc;
    }
    else return -2;

    df[0]+=dif1;
    df[1]+=dif2;
    if( (cabs(dif1)==0.0 && cabs(dif2)==0.0) || (cabs(dif1)<cabs(df[0])*reps && cabs(dif2)<cabs(df[1])*reps) ){
      *f*=I/(4.0*td->d);
      df[0]*=-1.0/(4.0*td->d);
      df[1]*=-( (signbit(r[1])==0)? 1.0 : -1.0)/(4.0*td->d);
      return 2*n+1;
    }
  }
  return -2;
}

int sr_qgf2(double complex *nf,double complex*df,double bn,double ay,double veps,double reps)
{
  double complex z,Fw,cep;
  double ep;
  
  z=bn*veps+I*ay/(2.0*veps);
  Fw=Faddeeva_w(z,reps);
  ep=exp(creal(z)*creal(z)-cimag(z)*cimag(z));
  cep=cos(bn*ay)+I*sin(bn*ay);
  
  *nf=2.0*I*ep*cimag(Fw)+2.0*cep;
  *df= -2.0*ep*creal(Fw)+2.0*cep;
  return 0;
}

int si_qgf2(double complex *nf,double complex*df,double cn,double ay,double veps,double reps)
{
  double cny,cne,aye,te,i_te,erp,erm,EXPLIMIT;

  EXPLIMIT=log(DBL_MAX);
  cny=cn*ay;
  cne=cn*veps;
  aye=ay/(2.0*veps);
  if(cny>EXPLIMIT){ // over flow
    if(-(cne+aye)*(cne+aye)+cny<-EXPLIMIT){ // under flow
      *nf=0.0;
      *df=0.0;
      return 0;
    }
    else return -1; // over flow
  }
  te=exp(cny);
  i_te=1.0/te;
  erp=erfc(cne+aye);
  erm=erfc(cne-aye);

  *nf= te*erp+i_te*erm;
  *df=-te*erp+i_te*erm;
  return 0;
}

int calc_cn(QPDT *td,QGEA *ta)
{
  long double pi2,tc,tp,tm,ti;
  
  pi2=8.0l*atanl(1.0l);
  tc=(long double)(td->d)/pi2;
  tp=tc*(long double)(td->kx+td->k);
  tm=tc*(long double)(td->kx-td->k);
  
  ta->nF[0]=(double)modfl(tp,&ti);  ta->nI[0]=(int)ti;
  ta->nF[1]=(double)modfl(tm,&ti);  ta->nI[1]=(int)ti;
  if(ta->nF[0]!=0.0 && ta->nF[1]!=0.0) return 0;
  else return -1;
}

double calc_bn2(int n,QGEA *ta)
{
  double dn,dnp,t2;

  dn=(double)(n-ta->nI[0])-ta->nF[0];
  dnp=dn*M_PI;
  t2=dnp+ta->c2;

  if(2.0*fabs(t2)/(fabs(dnp)+ta->c2) > 0.1) return -ta->c1*dn*t2;
  else {
    dn=(double)(n-ta->nI[1])-ta->nF[1];
    dnp=dn*M_PI;
    t2=dnp-ta->c2;
    return -ta->c1*dn*t2;
  }
}

int init_qgdt(QGDT *pa,double *r,QPDT *td)
{
  double tmp;

  pa->a =td->kx;
  pa->kx=td->k*r[0];
  pa->ky=td->k*r[1];
  pa->kd=td->k*td->d;
  pa->eikx[0]=cos(pa->kx)+I*sin(pa->kx);
  pa->eikx[1]=conj(pa->eikx[0]);
  tmp=-(td->k+td->kx)*td->d;
  pa->cel[0]=cos(tmp)+I*sin(tmp);
  tmp=-(td->k-td->kx)*td->d;
  pa->cel[1]=cos(tmp)+I*sin(tmp);

  pa->fc=0;
  return 0;
}

double complex tfu(double u,void *pa)
{
  QGDT *t=(QGDT *)pa;
  double complex sr,cc,ec1,ec2;
  double ekxu[2],ekdu;

  t->fc+=1;

  sr=csqrt(u*u-2.0*I*u);
  cc=ccos(fabs(t->ky)*sr);
  ekxu[0]=exp(t->kx*u);
  ekxu[1]=1.0/ekxu[0];
  ekdu=exp(t->kd*u);
  if(ekdu>DBL_MAX) return 0.0;
  if(exp((t->kd-fabs(t->kx))*u)>DBL_MAX) return 0.0;
  
  ec1=ekdu*t->cel[0]-1.0;
  ec2=ekdu*t->cel[1]-1.0;

  return (ekxu[1]*t->eikx[0]/ec1+ekxu[0]*t->eikx[1]/ec2)*cc/sr;
}

double complex tfx(double u,void *pa)
{
  QGDT *t=(QGDT *)pa;
  double complex sr,cc,ec1,ec2;
  double ekxu[2],ekdu;

  t->fc+=1;

  sr=csqrt(u*u-2.0*I*u);
  cc=ccos(fabs(t->ky)*sr);
  ekxu[0]=exp(t->kx*u);
  ekxu[1]=1.0/ekxu[0];
  ekdu=exp(t->kd*u);
  if(ekdu>DBL_MAX) return 0.0;
  if(exp((t->kd-fabs(t->kx))*u)>DBL_MAX) return 0.0;

  ec1=ekdu*t->cel[0]-1.0;
  ec2=ekdu*t->cel[1]-1.0;
  return (u-I)*(-ekxu[1]*t->eikx[0]/ec1+ekxu[0]*t->eikx[1]/ec2)*cc/sr;
}

double complex tfy(double u,void *pa)
{
  QGDT *t=(QGDT *)pa;
  double complex sr,cc,ec1,ec2;
  double ekxu[2],ekdu;

  t->fc+=1;

  sr=csqrt(u*u-2.0*I*u);
  cc=csin(t->ky*sr);
  ekxu[0]=exp(t->kx*u);
  ekxu[1]=1.0/ekxu[0];
  ekdu=exp(t->kd*u);
  if(ekdu>DBL_MAX) return 0.0;
  if(exp((t->kd-fabs(t->kx))*u)>DBL_MAX) return 0.0;

  ec1=ekdu*t->cel[0]-1.0;
  ec2=ekdu*t->cel[1]-1.0;
  return (ekxu[1]*t->eikx[0]/ec1+ekxu[0]*t->eikx[1]/ec2)*cc;
}

int qgf1_t1a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td)
{
  int s0_qgf1a(double *tf1,double *tf0,double kep,double  arg,double reps);
  int sl_qgf1a(double *af1,double *af0,double kep,double *apm,double reps); 

  double complex tep,epp,dif,dif0,dif1;
  double ra,xpr,ypr,kep,rp2,rm2,y2,tp1,tp0,af1[2],af0[2],apm[2],i_4e2,xp,xm,rp,rm,dr;
  int l,ef;

  ra=sqrt(r[0]*r[0]+r[1]*r[1]);
  xpr=r[0]/ra;
  ypr=r[1]/ra;
  
  kep=td->k*veps;
  i_4e2=1.0/(4.0*veps*veps);
  y2=r[1]*r[1];
  tep=cos(td->kx*td->d)+I*sin(td->kx*td->d);

  //  l=0;
  epp=1.0;
  rp2=r[0]*r[0]+y2;
  ef=s0_qgf1a(&tp1,&tp0,kep,rp2*i_4e2,reps);
  if(ef<0) return -1;
  else{
    *f=epp*tp1;
    df[0]=ra*epp*tp0;
    df[1]=0.0;
  }

  rp=ra;
  rm=ra;
  dr=xpr*td->d;
  xp=r[0];
  xm=r[0];
  // +-l
  for(l=1;l<QG1MAX;l++){
    xp+=td->d;
    xm-=td->d;
    rp+=dr;
    rm-=dr;
    epp*=tep;

    rp2=xp*xp+y2;
    rm2=xm*xm+y2;
    apm[0]=rp2*i_4e2;
    apm[1]=rm2*i_4e2;
    ef=sl_qgf1a(af1,af0,kep,apm,reps);
    if(ef<0) return -1;
    else {
      dif =epp*af1[0]+conj(epp)*af1[1];
      dif0=rp*epp*af0[0]+rm*conj(epp)*af0[1];
      dif1=(double)l*(epp*af0[0]-conj(epp)*af0[1]);

      *f+=dif;
      df[0]+=dif0;
      df[1]+=dif1;
      if( cabs(dif)<=cabs(*f)*reps && cabs(dif0)<=cabs(df[0])*reps && cabs(dif1)<=cabs(df[1])*reps ){
        *f/=4.0*M_PI;
        df[0]*= -1.0/(8.0*M_PI*veps*veps);
        df[1]*= td->d/(8.0*M_PI*veps*veps)*ypr;
        return 2*l+1;
      }
    }
  }
  return -1;
}

int s0_qgf1a(double *tf1,double *tf0,double kep,double arg,double reps)
{
  double i_fn,tfn,En1,En0,tkep2,kepn,expa,dif0,dif1;
  int n;

  expa=exp(-arg);
  i_fn=1.0;
  tkep2=kep*kep;
  kepn=1.0;
  if(arg>EXPINTLIMIT){
    *tf1=0.0;
    *tf0=0.0;
    return 0;
  }
  else {
    En1=expint(1,arg);
    En0=expint(0,arg);
  }
  *tf1=En1;
  *tf0=En0;
  for(n=1;n<QG1MAX;n++){
    tfn=1.0/(double)n;
    i_fn*=tfn;
    kepn*=tkep2;
    En0=En1;
    En1=tfn*(expa-arg*En0);
    dif0  =i_fn*kepn*En0;
    dif1  =i_fn*kepn*En1;
    *tf1+=dif1;
    *tf0+=dif0;
    if(dif1<0.0 || dif0<0.0) return n; // under flow
    else if( dif0<=*tf0*reps && dif1<=*tf1*reps ) return n; // convergent
  }
  return -1;
}

int sl_qgf1a(double *af1,double *af0,double kep,double *apm,double reps)
{
  double eap,eam,i_fn,tfn,tkep2,kepn,En1p,En1m,En0p,En0m,difp0,difm0,difp1,difm1;
  int n;
  
  eap=exp(-apm[0]);  eam=exp(-apm[1]);

  i_fn=1.0;
  tkep2=kep*kep;
  kepn=1.0;
  if(apm[0]>EXPINTLIMIT && apm[1]>EXPINTLIMIT){
    af1[0]=0.0;  af1[1]=0.0;
    af0[0]=0.0;  af0[1]=0.0;
    return 0;
  }
  else {
    En1p=expint(1,apm[0]);  En1m=expint(1,apm[1]);
    En0p=expint(0,apm[0]);  En0m=expint(0,apm[1]);
  }
  af1[0]=En1p;  af1[1]=En1m;
  af0[0]=En0p;  af0[1]=En0m;
  for(n=1;n<QG1MAX;n++){
    tfn=1.0/(double)n;
    i_fn*=tfn;
    kepn*=tkep2;
    En0p=En1p;
    En1p=tfn*(eap-apm[0]*En0p);
    En0m=En1m;
    En1m=tfn*(eam-apm[1]*En0m);
    difp0=i_fn*kepn*En0p;
    difm0=i_fn*kepn*En0m;
    difp1=i_fn*kepn*En1p;
    difm1=i_fn*kepn*En1m;
    af0[0]+=difp0;
    af0[1]+=difm0;
    af1[0]+=difp1;
    af1[1]+=difm1;
    if(difp0<0.0 || difm0<0.0 || difp1<0.0 || difm1<0.0) return n; // under flow
    else if( (difp0<=af0[0]*reps && difm0<=af0[1]*reps &&
              difp1<=af1[0]*reps && difm1<=af1[1]*reps) ) return n; // convergent
  }
  return -1;
}

int qgf1_t2a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td)
{
  double qgf1_inf(double u,void *pa);
  double qgf1_idf(double u,void *pa);
  double qgf1_mdf(double t,void *pa);

  QGEA tmpd;
  double complex epp,tep,difp,difm,it,dif,dif0,dif1;
  double ra,xpr,ypr,rl2,xp,xm,rp,rm,dr,y2,err,fpk2,i_ve2;
  int l;

  tmpd.c1=td->k*td->k;
  y2=r[1]*r[1];
  fpk2=4.0/tmpd.c1;
  i_ve2=1.0/(veps*veps);

  ra=sqrt(r[0]*r[0]+r[1]*r[1]);
  xpr=r[0]/ra;
  ypr=r[1]/ra;

  // l=0;
  epp=1.0;
  tep=cos(td->kx*td->d)+I*sin(td->kx*td->d);
  rl2=r[0]*r[0]+y2;
  tmpd.c2=0.25*rl2;
  
  *f=deintd(qgf1_inf,0.0,veps,&tmpd,reps,&err);  if(err<0.0) return -1;
  if(rl2>=fpk2){
    it=deintd(qgf1_idf,0.0,veps,&tmpd,reps,&err);  if(err<0.0) return -1;
  }
  else {
    it=2.0/rl2*deintid(qgf1_mdf,tmpd.c2*i_ve2,&tmpd,reps,&err);  if(err<0.0) return -1;
  }
  df[0]=ra*it;
  df[1]=0.0;

  rp=ra;
  rm=ra;
  dr=xpr*td->d;
  xp=r[0];
  xm=r[0];
  for(l=1;l<QG1MAX;l++){
    xp+=td->d;
    xm-=td->d;
    rp+=dr;
    rm-=dr;
    epp*=tep;
    // +l
    rl2=xp*xp+y2;
    tmpd.c2=0.25*rl2;
    dif=epp*deintd(qgf1_inf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    if(rl2>=fpk2){
      difp=epp*deintd(qgf1_idf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    }
    else {
      difp=epp*2.0/rl2*deintid(qgf1_mdf,tmpd.c2*i_ve2,&tmpd,reps,&err);  if(err<0.0) return -1;
    }
    // -l
    rl2=xm*xm+y2;
    tmpd.c2=0.25*rl2;
    dif+=conj(epp)*deintd(qgf1_inf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    if(rl2>=fpk2){
      difm=conj(epp)*deintd(qgf1_idf,0.0,veps,&tmpd,reps,&err); if(err<0.0) return -1;
    }
    else{
      difm=conj(epp)*2.0/rl2*deintid(qgf1_mdf,tmpd.c2*i_ve2,&tmpd,reps,&err);  if(err<0.0) return -1;
    }
    dif0=rp*difp+rm*difm;
    dif1=(double)l*(difp-difm);
    *f+=dif;
    df[0]+=dif0;
    df[1]+=dif1;
    if( cabs(dif)<=cabs(*f)*reps && cabs(dif0)<=cabs(df[0])*reps && cabs(dif1)<=cabs(df[1])*reps ){
      *f*=1.0/(2.0*M_PI);
      df[0]*= -1.0/(4.0*M_PI);
      df[1]*=td->d/(4.0*M_PI)*ypr;
      return 2*l+1;
    }
  }
  return -1;
}

double qgf1_mdf(double t,void *pa)
{
  QGEA *td=(QGEA *)pa;
  double te,i_t=1.0/t;;
  te=exp(-t+td->c1*td->c2*i_t);
  return te;
}

int qgf2_t1a(double complex *f,double complex *df,double *r,double veps,double reps,QPDT *td)
{
  int calc_cn(QPDT *td,QGEA *ta);
  double calc_bn2(int n,QGEA *ta);
  int sr_qgf2a(double complex *fn,double bn,double aye2,double aye1,double veps,double reps);
  int si_qgf2a(double *gn,double *hn,double cn,double ay,double aye1,double veps,double reps);

  QGEA ta;
  double complex eiax,fn,apb,eiby,dif,dif0,dif1;
  double kl,kh,k2,kd,tmp,anp,anm,bn,aye1,aye2,gn=0.0,hn=0.0,cn,ra,xpr,ypr,sgy;
  int nb,n;

  kl=0.8*td->k; // 
  kh=1.2*td->k; // 

  aye1=fabs(r[1])/(2.0*veps);
  aye2=aye1*aye1;
  ra=sqrt(r[0]*r[0]+r[1]*r[1]);
  xpr=r[0]/ra;
  ypr=r[1]/ra;
  sgy=((signbit(r[1])==0)? 1.0 : -1.0);
  
  k2=td->k*td->k;
  kd=2.0*M_PI/td->d;
  tmp=td->kx/kd;
  nb=(int)(tmp<0.0 ? tmp-0.5 : tmp+0.5);
  
  // init QGEA
  if(calc_cn(td,&ta)!=0) return -2;
  ta.c1=4.0*M_PI/(td->d*td->d);
  ta.c2=td->k*td->d;
  
  // n=nb;
  anp=kd*(double)nb-td->kx;
  anm=anp;
  if(fabs(anp)<kl || fabs(anp)>kh)   tmp=k2-anp*anp;
  else tmp=calc_bn2(nb,&ta);
  eiax=cos(anp*r[0])+I*sin(anp*r[0]);
  if(tmp>0.0){
    bn=sqrt(tmp);
    sr_qgf2a(&fn,bn,aye2,aye1,veps,reps);
    apb=anp/bn;
    eiby=cos(bn*fabs(r[1]))+I*sin(bn*fabs(r[1]));
    *f=eiax/bn*2.0*(I*cimag(fn)+eiby);
    df[0]=eiax*2.0*(I*xpr*apb*cimag(fn)-fabs(ypr)*creal(fn)+eiby*(xpr*apb+fabs(ypr)));
    df[1]=eiax*2.0*(I*ypr*apb*cimag(fn)+xpr*sgy*creal(fn)+eiby*(ypr*apb-xpr*sgy));
  }
  else if(tmp<0.0){
    cn=sqrt(-tmp);
    si_qgf2a(&gn,&hn,cn,fabs(r[1]),aye1,veps,reps);
    apb=anp/(I*cn);
    *f=eiax/(I*cn)*(gn + hn);
    df[0]=eiax*((xpr*apb+fabs(ypr))*gn+(xpr*apb-fabs(ypr))*hn);
    df[1]=eiax*((ypr*apb-xpr*sgy)*gn+(ypr*apb+xpr*sgy)*hn);
  }
  else return -2; // bn==0.0

  // nb+-n                                                                                                                                                                        
  for(n=1;n<QG2MAX;n++){
    anp+=kd;
    anm-=kd;
    // +
    if(fabs(anp)<kl || fabs(anp)>kh)   tmp=k2-anp*anp;
    else tmp=calc_bn2(nb+n,&ta);
    eiax=cos(anp*r[0])+I*sin(anp*r[0]);
    if(tmp>0.0){
      bn=sqrt(tmp);
      sr_qgf2a(&fn,bn,aye2,aye1,veps,reps);
      apb=anp/bn;
      eiby=cos(bn*fabs(r[1]))+I*sin(bn*fabs(r[1]));
      dif =eiax/bn*2.0*(I*cimag(fn)+eiby);
      dif0=eiax*2.0*(I*xpr*apb*cimag(fn)-fabs(ypr)*creal(fn)+eiby*(xpr*apb+fabs(ypr)));
      dif1=eiax*2.0*(I*ypr*apb*cimag(fn)+xpr*sgy*creal(fn)+eiby*(ypr*apb-xpr*sgy));
    }
    else if(tmp<0.0){
      cn=sqrt(-tmp);
      si_qgf2a(&gn,&hn,cn,fabs(r[1]),aye1,veps,reps);
      apb=anp/(I*cn);
      dif =eiax/(I*cn)*(gn+hn);
      dif0=eiax*((xpr*apb+fabs(ypr))*gn+(xpr*apb-fabs(ypr))*hn);
      dif1=eiax*((ypr*apb-xpr*sgy)*gn+(ypr*apb+xpr*sgy)*hn);
    }
    else return -2;
    
    // -
    if(fabs(anm)<kl || fabs(anm)>kh)   tmp=k2-anm*anm;
    else tmp=calc_bn2(nb-n,&ta);
    eiax=cos(anm*r[0])+I*sin(anm*r[0]);
    if(tmp>0.0){
      bn=sqrt(tmp);
      sr_qgf2a(&fn,bn,aye2,aye1,veps,reps);
      apb=anm/bn;
      eiby=cos(bn*fabs(r[1]))+I*sin(bn*fabs(r[1]));
      dif +=eiax/bn*2.0*(I*cimag(fn)+eiby);
      dif0+=eiax*2.0*(I*xpr*apb*cimag(fn)-fabs(ypr)*creal(fn)+eiby*(xpr*apb+fabs(ypr)));
      dif1+=eiax*2.0*(I*ypr*apb*cimag(fn)+xpr*sgy*creal(fn)+eiby*(ypr*apb-xpr*sgy));
    }
    else if(tmp<0.0){
      cn=sqrt(-tmp);
      si_qgf2a(&gn,&hn,cn,fabs(r[1]),aye1,veps,reps);
      apb=anm/(I*cn);
      dif +=eiax/(I*cn)*(gn+hn);
      dif0+=eiax*((xpr*apb+fabs(ypr))*gn+(xpr*apb-fabs(ypr))*hn);
      dif1+=eiax*((ypr*apb-xpr*sgy)*gn+(ypr*apb+xpr*sgy)*hn);
    }
    else return -2;
    
    *f+=dif;
    df[0]+=dif0;
    df[1]+=dif1;
    if(  cabs(dif)<=cabs(*f)*reps && cabs(dif0)<=cabs(df[0])*reps && cabs(dif1)<=cabs(df[1])*reps ){
      *f*=I/(4.0*td->d);
      df[0]*=-1.0/(4.0*td->d);
      df[1]*= 1.0/(4.0*td->d);
      return 2*n+1;
    }
  }
  return -2;
}

int sr_qgf2a(double complex *fn,double bn,double aye2,double aye1,double veps,double reps)
{
  *fn=exp(bn*bn*veps*veps-aye2)*Faddeeva_w(bn*veps+I*aye1,reps); 
  return 0;
}

int si_qgf2a(double *gn,double *hn,double cn,double ay,double aye,double veps,double reps)
{
  double EXPLIMIT,cny,cne,te,i_te,erp,erm;

  EXPLIMIT=log(DBL_MAX);
  cny=cn*ay;
  cne=cn*veps;
  if(cny>EXPLIMIT){ // over flow
    if(-(cne+aye)*(cne+aye)+cny<-EXPLIMIT){ // under flow
      *gn=0.0;
      *hn=0.0;
      return 0;
    }
    else return -1; // over flow
  }
  te=exp(cny);
  i_te=1.0/te;
  erp=erfc(cne+aye);
  erm=erfc(cne-aye);

  *gn=i_te*erm;
  *hn=te*erp;
  return 0;
}

