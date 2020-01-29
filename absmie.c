/* absmie.c

This is a program for Mie scattering in absorbing medium.

-------------------------------------------------------------------
I Wayan Sudiarta, 
Department of Physics and Atmospheric Science
Dalhousie University, Halifax, NS, Canada B3H 3J5
sudiarta@dal.ca

-------------------------------------------------------------------
References:

I.W. Sudiarta and P. Chylek, "Mie-scattering formalism for spherical particles embedded in an absorbing medium", J. Opt. Soc. Am. A 18 (6): 1275-1278 (2001).

I.W. Sudiarta, " Effective medium approximation for light scattering of heterogeneous particles", PhD Thesis, Dalhousie University, 2003.

---------------------- --------------------------------------------
  Original program taken from Bohren and Huffman (1983), Appendix A
  Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
  in order to compute <cos(theta)>
  
  This code was translatted to C by P. J. Flatau Feb 1998. The C
  version uses "Numerical Recipes" public domain code for complex
  arithmetics "complex.c" and "nrutil.c" (http://www.nr.com).
-------------------------------------------------------------------

-- November 15, 1999 --

  And modified by I Wayan Sudiarta 
  for calculating a scattering in absorbing medium
 
-- 27 August 2001 --
  Adding Phase Function and Asymmetry

 
*/

#include <math.h>
#include <stdio.h>
#include "complex.h"
#include "nrutil.h"

#define nmxx 5000
#define nangel 500
#define CXONE Complex(1.0, 0.0)  /* Complex 1 */
#define PI 4.E0*atan(1.E0)
#define CXI Complex(0.0,1.0E0)  /* Complex i */

void mie(float a, fcomplex k, fcomplex k1, float  *qext,
   float *qsca,float *qabs, float ph[nangel], float *g,float *ii);

int main(void)
{
  int i;
  float gg, dangle, qext, qsca, qabs, ii;
  float x; 		/* size parameter 2 Pi Radius of sphere /lambda */
  fcomplex medref;      /* complex medium */
  fcomplex sphereref;   /* complex sphere */
  float phase[nangel];  /* phase functions */
  float tr, ti;

  FILE *infile, *outfile;

  dangle = 180/(float)(nangel-1);

  // Read Input File  
  infile = fopen("mie.dat","r");
  if(infile==NULL){printf("please provide mie.dat\n"); exit(1);}
  fscanf(infile,"%f",&x);
  fscanf(infile,"%f  %f",&tr,&ti);
  sphereref = Complex(tr,ti);
  fscanf(infile,"%f  %f",&tr,&ti);
  medref = Complex(tr, ti);
  fclose(infile);

  // Compute   
  mie(x, medref, sphereref, &qext, &qsca, &qabs, phase, &gg, &ii);

  // Print Output
  printf(" Mie Scattering in Absorbing Medium \n");
  printf(" size parameter = %f\n",x);
  printf(" msphere        = %f  %f \n",sphereref.r,sphereref.i);
  printf(" mmed           = %f  %f \n",medref.r,medref.i);

  printf("\n Qext, Qsca, Qabs, g, Ii, %f  %f  %f  %f  %f\n\n",qext,qsca,qabs,gg,ii);
  printf("\n Angle......phase \n");
  for(i=0;i<nangel;i++){
      printf("%f %f \n",i*dangle, phase[i]);
  }

  // Save Output
  outfile = fopen("mie.out","w");
  // Print Output
  fprintf(outfile," Mie Scattering in Absorbing Medium \n");
  fprintf(outfile," size parameter = %f\n",x);
  fprintf(outfile," msphere        = %f  %f \n", sphereref.r, sphereref.i);
  fprintf(outfile," mmed           = %f  %f \n",medref.r,medref.i);

  fprintf(outfile,"\n Qext, Qsca, Qabs, g, Ii, %f %f %f %f %f\n\n",qext,qsca,qabs,gg,ii);
  fprintf(outfile,"\n Angle......phase \n");
  for(i=0;i<nangel;i++){
      fprintf(outfile,"%f %f \n",i*dangle, phase[i]);
  }
  fclose(outfile);

  return 0;
}

void  mie(float a, fcomplex m, fcomplex m1, float  *qext, float *qsca,
float *qabs, float ph[nangel], float *g, float *ii)
{
      unsigned int j, n, nmx, nn, nstop;
      float pii, t, theta, xstop, ymod;
      fcomplex cxan, cxan1, cxbn, cxbn1, cxxi, cxxi0, cxy, cxxi1;
      fcomplex x, y, cxref, cxtemp, cxcn, cxdn;
      fcomplex in, chi, chi0, chi1;
      fcomplex temp1, temp2, temp3, temp4, cxan2, cxbn2;
      fcomplex qexttemp,qabstemp, qscatemp;
      fcomplex psi, psi0, psi1, psitemp, xitemp;
      fcomplex psi2, psi20, psi21;
      fcomplex s1[nangel],s2[nangel];
      float mu[nangel], pi[nangel], pi0[nangel], pi1[nangel],tau[nangel];
      float dang, fn, phtemp, fn1, fn2, phtemp1, phtemp2, gtemp;      

      double  rn, A; 


/* .. Local Arrays ..*/
      fcomplex cxd[nmxx];

      x = RCmul(a,m);
      cxref = Cdiv(m1,m);
      y = RCmul(a,m1);

/* . . .  Calculating A  . . . */
      if(x.i != 0.0){
        A = m.r*(exp(2.0* x.i)/(2.0* x.i) 
	    + (1.0-exp(2.0*x.i))/((2.0*x.i)*(2.0* x.i))); 
      }
      else{
        A = 0.5*m.r; 
      }
      
      *ii = A;

/* Series expansion terminated after NSTOP terms */
      xstop = x.r + 4.E0*pow(x.r,0.3333) + 2.0;
      nstop = xstop;
      ymod = Cabs(y);
      nmx = FMAX(xstop,ymod) + 15;

      if (nmx>nmxx) {
        printf(" x, nmx, nmxx, cxref %f %i %i  \n ", x.r, nmx, nmxx);
        printf(" xstop nstop ymod %f %i %f \n", xstop, nstop, ymod); 
        printf(" Error: NMX > NMXX= %i \n", nmxx);
        return;
      }
      
    /* calculate the initial iteration for calculating pi and tau */

      dang = PI/(float)(nangel-1);
      for (j = 0; j < nangel; j++) 
      {
        theta = (float)(j)*dang;
        mu[j] = cos(theta);
      }
        
     for ( j = 0; j < nangel; j++) 
     {
        pi0[j] = 0.E0;
        pi1[j] = 1.E0;
	s1[j] = Complex(0.0,0.0);
	s2[j] = Complex(0.0,0.0);
     }


/* Logarithmic derivative D(J) calculated by downward recurrence
    beginning with initial value (0.,0.) at J=NMX */

      cxd[nmx] = Complex(0.E0,0.E0);
      nn = nmx - 1;

      for (n = 1; n<= nn; n++) {
        rn = nmx - n + 1;
/*        cxd[nmx-n] = (rn/y) - (1.E0/(cxd[nmx-n+1]+rn/y)) */
        cxtemp=Cadd(cxd[nmx-n+1],Cdiv(Complex(rn,0.0),y));
        cxtemp=Cdiv(CXONE,cxtemp);
        cxd[nmx-n]=Csub(Cdiv(Complex(rn,0.0),y),cxtemp);
      }

/* Riccati-Bessel functions with complex argument X
    calculated by upward recurrence - modified by Wayan Nov 3*/

      psi0 = Ccos(x);
      psi1 = Csin(x);
      chi0 = RCmul(-1.0,Csin(x));
      chi1 = Ccos(x);

      psi20 = Ccos(y);
      psi21 = Csin(y);

      /* apsi - i* chi */
      cxxi0 = Complex(psi0.r + chi0.i, psi0.i - chi0.r);
      cxxi1 = Complex(psi1.r + chi1.i, psi1.i - chi1.r);
      
      qexttemp = Complex(0.E0,0.E0);
      qabstemp = Complex(0.E0,0.E0);	 
      qscatemp = Complex(0.E0,0.E0);
      phtemp = 0.0E0;

      in = Complex(0.0, 1.0E0);

/* == summations == */
      
for ( n = 1; n <= nstop; n++) {  
  rn = n;
  psi = Cadd(Cdiv(RCmul(2.E0*rn-1.E0,psi1),x), RCmul(-1.E0,psi0));
  psi2 = Cadd(Cdiv(RCmul(2.E0*rn-1.E0,psi21),y),RCmul(-1.E0,psi20)); 
  chi = Cadd(Cdiv(RCmul(2.E0*rn-1.E0,chi1),x), RCmul(-1.E0,chi0));
  cxxi = Complex(psi.r + chi.i, psi.i - chi.r);

/* Compute AN and BN:*/
/*        cxan = (cxd(n)/cxref+rn/x)*psi - psi1; */

  cxan=Cdiv(cxd[n],cxref);
  cxan=Cadd(cxan,Cdiv(Complex(rn,0.0),x));
  cxan=Cmul(cxan,psi);
  cxan=Csub(cxan,psi1);

/*        cxan = cxan/((cxd(n)/cxref+rn/x)*cxxi-cxxi1); */
  cxtemp=Cdiv(cxd[n],cxref);
  cxtemp=Cadd(cxtemp,Cdiv(Complex(rn,0.0),x));
  cxtemp=Cmul(cxtemp,cxxi);
  cxtemp=Csub(cxtemp,cxxi1);
  cxan=Cdiv(cxan,cxtemp);
	
 cxdn=Cdiv(Complex(0.0,-1.E0),Cmul(psi2,cxtemp)); 

/*        cxbn = (cxref*cxd(n)+rn/x)*psi - psi1; */
  cxbn=Cmul(cxref,cxd[n]);
  cxbn=Cadd(cxbn,Cdiv(Complex(rn,0.0),x));
  cxbn=Cmul(cxbn,psi);
  cxbn=Csub(cxbn,psi1);

/*        cxbn = cxbn/((cxref*cxd(n)+rn/x)*cxxi-cxxi1); */
  cxtemp=Cmul(cxref,cxd[n]);
  cxtemp=Cadd(cxtemp,Cdiv(Complex(rn,0.0),x));
  cxtemp=Cmul(cxtemp,cxxi);
  cxtemp=Csub(cxtemp,cxxi1);
  cxbn=Cdiv(cxbn,cxtemp);

  cxcn=Cdiv(Cmul(Complex(0.0,-1.E0),cxref),Cmul(cxtemp,psi2)); 


  psitemp = Cadd(psi1,RCmul(-1.E0*rn,Cdiv(psi,x)));
  xitemp =Cadd(cxxi1,RCmul(-1.E0*rn,Cdiv(cxxi,x)));
 
  temp1 = Cmul(xitemp,Conjg(cxxi)); /* xi'xi* */
  temp2 = Cmul(xitemp,Conjg(psi));  /* xi'psi* */
  temp3 = Cmul(psitemp,Conjg(cxxi)); /* psi'xi* */
  temp4 = Cmul(psitemp,Conjg(psi)); /* psi'psi* */

  /* calculating the extinction */

  cxtemp = Cmul(Complex(0.E0,1.E0),temp4);
  cxtemp = Cadd(cxtemp, Cmul(Complex(0.E0,-1.E0),Conjg(temp4)));
  cxtemp = Cadd(cxtemp, Cmul(cxbn, Cmul(Complex(0.E0,1.E0),Conjg(temp3))));
  cxtemp = Cadd(cxtemp, Cmul(Conjg(cxbn), Cmul(Complex(0.E0,1.E0),Conjg(temp2))));
  cxtemp = Cadd(cxtemp, Cmul(cxan, Cmul(Complex(0.E0,-1.E0),temp2)));
  cxtemp = Cadd(cxtemp, Cmul(Conjg(cxan), Cmul(Complex(0.E0,-1.E0),temp3)));

  qexttemp = Cadd(qexttemp, RCmul(2*rn+1.0, cxtemp));

  /* scattering crosssection */

  cxtemp = Cmul(RCmul(Cabs(cxan)*Cabs(cxan),Complex(0.E0,-1.E0)),temp1);
  cxtemp = Cadd(cxtemp, Cmul(RCmul(Cabs(cxbn)*Cabs(cxbn),Complex(0.E0,1.E0)),Conjg(temp1)));

  qscatemp = Cadd(qscatemp, RCmul(2*rn+1, cxtemp));

  /* absorption  */
  psitemp =Cadd(psi21,RCmul(-1.E0*rn,Cdiv(psi2,y)));
  temp1 = Cmul(Conjg(psi2), psitemp); /* psi2'psi2* */
  cxtemp = Cmul(RCmul(Cabs(cxdn)*Cabs(cxdn),Complex(0.E0,1.E0)),temp1);
  cxtemp = Cadd(cxtemp, Cmul(RCmul(Cabs(cxcn)*Cabs(cxcn),Complex(0.E0,-1.E0)),Conjg(temp1)));

  qabstemp = Cadd(qabstemp, RCmul(2*rn+1, cxtemp));

  fn = (2.E0*rn+1.E0)/(rn*(rn+1.E0));
  
  /* calculate phase temp */

  phtemp1 = Cabs(cxan);
  phtemp2 = Cabs(cxbn);
  phtemp = phtemp + (2.0*rn+1.0)*(phtemp1*phtemp1 + phtemp2*phtemp2);

  /* calculate Asymmetry */
  if(n==1){
   cxan2 = cxan;
   cxbn2 = cxbn;
   gtemp = 0.0;
  }
  else {
    fn1 = (rn-1.0)*(rn+1.0E0)/rn;
    fn2 = (2.E0*(rn-1.0)+1.E0)/((rn-1.0)*rn);
    temp1 = Cadd(Cmul(cxan2,Conjg(cxan)),Cmul(cxbn2,Conjg(cxbn)));
    temp2 = Cmul(cxan2,Conjg(cxbn2));
    gtemp = gtemp + fn1*temp1.r + fn2*temp2.r;
    cxan2 = cxan;
    cxbn2 = cxbn;
  }       

/* Calculate Phase Function */

for(j=0;j<nangel;j++)
{

/* updating  pi, tau */

    pi[j] = pi1[j];
    tau[j] = rn*mu[j]*pi[j] - (rn+1.E0)*pi0[j];

/* calculation of s1 and s2 */

    temp1 = Cadd(RCmul(pi[j],cxan),RCmul(tau[j],cxbn));
    s1[j] = Cadd(s1[j],RCmul(fn,temp1));
    
    temp2 = Cadd(RCmul(tau[j],cxan),RCmul(pi[j],cxbn));
    s2[j] = Cadd(s2[j],RCmul(fn,temp2));   
   
}

  /*  compute pi_n+1 from PI = pi_n , PI0 = pi_n-1 */
            
    for ( j = 0; j< nangel; j++) {
       pi1[j] = ((2.*rn+1.)*mu[j]*pi[j]-(rn+1.)*pi0[j])/rn;
       pi0[j] = pi[j];
    }
        
    in = Cmul(in,Complex(0.0, 1E0));

    psi0 = psi1;
    psi1 = psi;
    psi20 = psi21;
    psi21 = psi2;
    chi0 = chi1;
    chi1 = chi;
    cxxi1 = Complex(psi1.r + chi1.i, psi1.i - chi1.r);

} /*end of big for */

  /*  Have summed sufficient terms Now compute *qsca,*qext,*qabs */

  qexttemp = Cmul(Conjg(m),qexttemp);
  qscatemp = Cmul(Conjg(m),qscatemp);
  qabstemp = Cmul(Conjg(m1),qabstemp);
    
  *qsca = (1.E0/(A*Cabs(x)*Cabs(x)))* qscatemp.r;
  *qext = (1.E0/(A*Cabs(x)*Cabs(x)))* qexttemp.r;
  *qabs = (1.E0/(A*Cabs(y)*Cabs(y)))* qabstemp.r;
  
  /* Phase function */
  for(j=0;j<nangel;j++){
     phtemp1 = Cabs(s1[j]);
     phtemp2 = Cabs(s2[j]);
     ph[j] = (phtemp1*phtemp1 + phtemp2*phtemp2)/phtemp;
  }
  
  *g = 2.0*gtemp/phtemp;

}



