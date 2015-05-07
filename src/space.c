/*
##
## Copyright (c) 2011, Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission. 
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "zufall.h"
//#include "matrix.h"
//#include "mxs.h"
#include "baymri.h"

/******************************/
/* Global Variables			  */
/******************************/
//int setprint = 1; // Start proc()... printed for all procedures
// int setprint = 0;

void choldc(double **a, double **res, int n, int bw, int start)
{ 
 int i,j,k;
 double sum;
 for (i=start;i<=n;i++)
   {
     for (j=0;j<=bw;j++) 
       {
	 for (sum=a[i][j],k=i-1;(i-k<=bw)&&(k>=1) ;k--) 
	   {
	     sum -= res[k][i-k]*res[k][j];
	   }
	 if (0 == j) 
	   {
	     res[i][0]=sqrt(sum);
	   }
	 else res[i][j]=sum/res[i][0];
       }
   }
}

double logit(double x)
{
  return(log(x/(1-x)));
}

double power(double x, int n)
{
  if (n==0){return 1;}
  if (n==1){return x;}
  return (x*power(x,n-1));
}

int ix(int x,int y,int z,int X,int Y, int Z)
{
  return((z-1)*Y*X+(y-1)*X+x-1);
}

int indx(int x,int y,int z, int t, int X, int Y, int Z, int T) 
{
  return((t-1)*X*Y*Z+(z-1)*X*Y+(y-1)*X+x-1);
}

double time1_last_time=0.0;
int time1_last_result=0;
int time1(double t,double* CPtime)
{
  int start=0;
  double zeit=0.0;
  if (time1_last_time==t)
    {
      start=time1_last_result;
    }
  else
    {
      if (time1_last_time<t){
	start=time1_last_result;
	zeit=time1_last_time;
      }
      
      while (zeit<t)
	{
	  start++;
	  zeit=CPtime[start];
	}
      time1_last_time=t;
      time1_last_result=start;
    }
  return start;
}

double space_log_fc_gamma(double gamma, double tau_epsilon, double tau_gamma, double* conc, double* time, double t0, double kep, double vp, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings)
 {
   double p = 0.0;
   p -= 0.5*tau_gamma*gamma*gamma;
   for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,time[t]+t0,aifsettings)-exp(gamma)*convterm(kep,time[t]+t0,aifsettings),2);
     }
   return p;
 }

double space_log_fc_gamma2(double gamma, double* ktrans, double tau_epsilon, double* tau_gamma, double* tau_gamma2, double* conc, double* time, double t0, double kep, double vp, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings, int* img_mask)
 {
   double p = 0.0;
   if (x!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma[ix(x-1,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x-1,y,z,X,Y,Z)]),2.0);}}
   if (x!=X){if(img_mask[ix(x+1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma[ix(x,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x+1,y,z,X,Y,Z)]),2.0);}}
   if (y!=1){if(img_mask[ix(x,y-1,z,X,Y,Z)]){p -= 0.5*tau_gamma2[ix(x,y-1,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y-1,z,X,Y,Z)]),2.0);}}
   if (y!=Y){if(img_mask[ix(x,y+1,z,X,Y,Z)]){p -= 0.5*tau_gamma2[ix(x,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y+1,z,X,Y,Z)]),2.0);}}
   for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,time[t]+t0,aifsettings)-exp(gamma)*convterm(kep,time[t]+t0,aifsettings),2);
     }
   return p;
 }

double space_log_fc_gamma3(double gamma, double* ktrans, double tau_epsilon, double* tau_gamma, double* tau_gamma2, double* tau_gamma3, double* conc, double* time, double t0, double kep, double vp, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings, int* img_mask)
 {
   double p = 0.0;
   if (x!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma[ix(x-1,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x-1,y,z,X,Y,Z)]),2.0);}}
   if (x!=X){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma[ix(x,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x+1,y,z,X,Y,Z)]),2.0);}}
   if (y!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma2[ix(x,y-1,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y-1,z,X,Y,Z)]),2.0);}}
   if (y!=Y){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma2[ix(x,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y+1,z,X,Y,Z)]),2.0);}}
   if (z!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma3[ix(x,y,z-1,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y,z-1,X,Y,Z)]),2.0);}}
   if (z!=Z){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma3[ix(x,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y,z+1,X,Y,Z)]),2.0);}}
   for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,time[t]+t0,aifsettings)-exp(gamma)*convterm(kep,time[t]+t0,aifsettings),2);
     }
   return p;
 }


double space_log_fc_theta(double theta, double tau_epsilon, double tau_theta, double* conc, double* time, double t0, double ktrans, double vp, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings)
 {
   double p = 0.0;
   p -= 0.5*tau_theta*theta*theta;
   for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,time[t]+t0,aifsettings)-ktrans*convterm(exp(theta),time[t]+t0,aifsettings),2);
     }
   return p;
 }

double log_fc_eta(double eta, double tau_epsilon, double tau_eta, double* conc, double* time, double t0, double kep, double ktrans, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings)
 {
   double p = 0.0;
   p -= 0.5*tau_eta*eta*eta;
   double a,b;
   for (int t=1; t<=T; t++)
     {
       a = aif(time[t]+t0, aifsettings);
       b=a*(conc[indx(x,y,z,t,X,Y,Z,T)]-ktrans*convterm(kep,time[t]+t0,aifsettings));
       p -= 0.5*tau_epsilon*(a*exp(2*eta/(1+eta))-2*b*exp(eta/(1+eta)));
     }

   return p;
 }

double space_log_fc_eta3(double eta, double tau_epsilon, double tau_eta, double* conc, double* time, double t0, double kep, double ktrans, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings, double a_vp, double b_vp)
 {
   double p = 0.0;
   double a,b;
   for (int t=1; t<=T; t++)
     {
       a = aif(time[t]+t0,aifsettings);
       b=conc[indx(x,y,z,t,X,Y,Z,T)]-ktrans*convterm(kep,time[t]+t0,aifsettings);
       p -= 0.5*tau_epsilon*(a*eta-b)*(a*eta-b);
     }
   p += (a_vp-1)*log(eta)+(b_vp-1)*log(1-eta);
   return p;
 }


double log_fc_eta2(double eta, double* vp, double tau_epsilon, double* tau_eta, double* tau_eta2, double* conc, double* time, double t0, double kep, double ktrans, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings)
 {
   double p = 0.0;
   if (x!=1){p -= 0.5*tau_eta[ix(x-1,y,z,X,Y,Z)]*pow(logit(vp[ix(x-1,y,z,X,Y,Z)]),2.0);}
   if (x!=X){p -= 0.5*tau_eta[ix(x,y,z,X,Y,Z)]*pow(logit(vp[ix(x+1,y,z,X,Y,Z)]),2.0);}
   if (y!=1){p -= 0.5*tau_eta2[ix(x,y-1,z,X,Y,Z)]*pow(logit(vp[ix(x,y-1,z,X,Y,Z)]),2.0);}
   if (y!=Y){p -= 0.5*tau_eta2[ix(x,y,z,X,Y,Z)]*pow(logit(vp[ix(x,y+1,z,X,Y,Z)]),2.0);}
   double a,b;
   for (int t=1; t<=T; t++)
     {
       a = aif(time[t]+t0,aifsettings);
       b=a*(conc[indx(x,y,z,t,X,Y,Z,T)]-ktrans*convterm(kep,time[t]+t0,aifsettings));
       p -= 0.5*tau_epsilon*exp((2.0*log(a)-log(2.0*b))*eta/(1+eta));
   }

   return p;
 }

double space_log_fc_theta_tau(double tau_theta, double theta, double tau_epsilon, double ktrans, double* conc, double* time, double t0, double vp, int x, int y, int z, int X, int Y, int Z, int T, double aa, double bb, double* aifsettings)
 {
   double p = 0.0;
   p += log(tau_theta-sqrt(2*PI))+0.5*tau_theta*theta*theta;
   p += (aa*log(tau_theta))-bb*tau_theta;

       for (int t=1; t<=T; t++)
	 {
	   p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,time[t]+t0,aifsettings)-ktrans*convterm(exp(theta),time[t]+t0,aifsettings),2);
	 }
   return p;
 }

double space_log_fc_gamma_tau(double tau_gamma,double gamma, double tau_epsilon, double kep, double* conc, double* time, double t0, double vp, int x, int y, int z, int X, int Y, int Z, int T, double aa, double bb, double* aifsettings)
 {
   double p = 0.0;
   p += log(tau_gamma-sqrt(2*PI))+0.5*tau_gamma*gamma*gamma;
   p += (aa*log(tau_gamma))-bb*tau_gamma;
   for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,time[t]+t0,aifsettings)-exp(gamma)*convterm(kep,time[t]+t0,aifsettings),2);
     }
   return p;
 }


double space_update_gamma2(double gamma, double* ktrans, double kep, double vp, double* tau_gamma, double* tau_gamma2, double* tau_gamma3,double tau_epsilon, double* conc, double* time, double t0, int x, int y, int z, int X, int Y, int Z, int T, double sigma, double* aifsettings, int spatial, int* img_mask)
{

//	if(setprint==1){printf("Start space_update_gamma2() \n");}

  double gamma_new = normal(gamma,sigma);
    
  double logalpha = 0.0;

    if (spatial==3)
      {
	logalpha += space_log_fc_gamma3(gamma_new,ktrans,tau_epsilon,tau_gamma,tau_gamma2,tau_gamma3,conc,time,t0,kep,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
	logalpha -= space_log_fc_gamma3(gamma,ktrans,tau_epsilon,tau_gamma,tau_gamma2,tau_gamma3,conc,time,t0,kep,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
      }
    if (spatial==1)
      {
	logalpha += space_log_fc_gamma2(gamma_new,ktrans,tau_epsilon,tau_gamma,tau_gamma2,conc,time,t0,kep,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
	logalpha -= space_log_fc_gamma2(gamma,ktrans,tau_epsilon,tau_gamma,tau_gamma2,conc,time,t0,kep,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
      }
  if (spatial==0)
    {
      logalpha +=  space_log_fc_gamma(gamma_new,tau_epsilon,tau_gamma[ix(x,y,z,X,Y,Z)],conc,time,t0,kep,vp,x,y,z,X,Y,Z,T,aifsettings);
      logalpha -= space_log_fc_gamma(gamma,tau_epsilon,tau_gamma[ix(x,y,z,X,Y,Z)],conc,time,t0,kep,vp,x,y,z,X,Y,Z,T,aifsettings);
    }
  if (exp(logalpha)>nulleins())
    {
      return gamma_new;
    }
  else
    {
      return gamma;
    }
}


double space_update_eta2(double eta, double kep, double ktrans, double* tau_eta, double* tau_eta2, double tau_epsilon, double* conc, double* time, double t0, int x, int y, int z, int X, int Y, int Z, int T, double sigma, double* aifsettings)
{
//	if(setprint==1){Rprintf("Start space_update_eta2() \n");}

  double eta_new = normal(eta,sigma);
  double logalpha = 0.0;

 logalpha +=  log_fc_eta(eta_new,tau_epsilon,tau_eta[ix(x,y,z,X,Y,Z)],conc,time,t0,kep,ktrans,x,y,z,X,Y,Z,T,aifsettings);
 logalpha -= log_fc_eta(eta,tau_epsilon,tau_eta[ix(x,y,z,X,Y,Z)],conc,time,t0,kep,ktrans,x,y,z,X,Y,Z,T,aifsettings);

  if (exp(logalpha)>nulleins())
    {
      return eta_new;
    }
  else
    {
      return eta;
    }
}

double space_update_eta3(double vp, double kep, double ktrans, double* tau_eta, double* tau_eta2, double tau_epsilon, double* conc, double* time, double t0, int x, int y, int z, int X, int Y, int Z, int T, double sigma, double* aifsettings,double a_vp, double b_vp)
{
//	if(setprint==1){Rprintf("Start space_update_eta3() \n");}

  double eta_new = -1;
  while (eta_new>1 || eta_new<0)
    {
      eta_new=normal(vp,sigma);
    }
  double logalpha = 0.0;

 logalpha +=  space_log_fc_eta3(eta_new,tau_epsilon,tau_eta[ix(x,y,z,X,Y,Z)],conc,time,t0,kep,ktrans,x,y,z,X,Y,Z,T,aifsettings,a_vp,b_vp);
 logalpha -= space_log_fc_eta3(vp,tau_epsilon,tau_eta[ix(x,y,z,X,Y,Z)],conc,time,t0,kep,ktrans,x,y,z,X,Y,Z,T,aifsettings,a_vp,b_vp);
  if (exp(logalpha)>nulleins())
    {
      return eta_new;
    }
  else
    {
      return vp;
    }
}

double update_tau(double tau, double a, double b, double theta)
{
//	if(setprint==1){Rprintf("Start update_tau() \n");}
  return(RNDGAM(a+1,b+0.5*pow(theta,2.0)));
}


double space_log_fc_theta2(double theta, double* kep, double tau_epsilon, double* tau_theta, double* tau_theta2, double* conc, double* time, double t0, double ktrans, double vp, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings, int* img_mask)
 {
   double p = 0.0;
   if (x!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_theta[ix(x-1,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x-1,y,z,X,Y,Z)]),2.0);}}
   if (x!=X){if(img_mask[ix(x+1,y,z,X,Y,Z)]){p -= 0.5*tau_theta[ix(x,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x+1,y,z,X,Y,Z)]),2.0);}}
   if (y!=1){if(img_mask[ix(x,y-1,z,X,Y,Z)]){p -= 0.5*tau_theta2[ix(x,y-1,z,X,Y,Z)]*pow(theta-log(kep[ix(x,y-1,z,X,Y,Z)]),2.0);}}
   if (y!=Y){if(img_mask[ix(x,y+1,z,X,Y,Z)]){p -= 0.5*tau_theta2[ix(x,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x,y+1,z,X,Y,Z)]),2.0);}}
    for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,time[t]+t0,aifsettings)-ktrans*convterm(exp(theta),time[t]+t0,aifsettings),2);
     }
   return p;
 }

double space_log_fc_theta3(double theta, double* kep, double tau_epsilon, double* tau_theta, double* tau_theta2, double* tau_theta3,double* conc, double* ttime, double t0, double ktrans, double vp, int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings, int* img_mask)
 {
   double p = 0.0;
   if (x!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_theta[ix(x-1,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x-1,y,z,X,Y,Z)]),2.0);}}
   if (x!=X){if(img_mask[ix(x+1,y,z,X,Y,Z)]){p -= 0.5*tau_theta[ix(x,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x+1,y,z,X,Y,Z)]),2.0);}}
   if (y!=1){if(img_mask[ix(x,y-1,z,X,Y,Z)]){p -= 0.5*tau_theta2[ix(x,y-1,z,X,Y,Z)]*pow(theta-log(kep[ix(x,y-1,z,X,Y,Z)]),2.0);}}
   if (y!=Y){if(img_mask[ix(x,y+1,z,X,Y,Z)]){p -= 0.5*tau_theta2[ix(x,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x,y+1,z,X,Y,Z)]),2.0);}}
   if (z!=1){if(img_mask[ix(x,y,z-1,X,Y,Z)]){p -= 0.5*tau_theta3[ix(x,y,z-1,X,Y,Z)]*pow(theta-log(kep[ix(x,y,z-1,X,Y,Z)]),2.0);}}
   if (z!=Z){if(img_mask[ix(x,y,z+1,X,Y,Z)]){p -= 0.5*tau_theta3[ix(x,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x,y,z+1,X,Y,Z)]),2.0);}}
    for (int t=1; t<=T; t++)
     {
       p -= 0.5*tau_epsilon*power((conc[indx(x,y,z,t,X,Y,Z,T)])-extraterm(vp,ttime[t]+t0,aifsettings)-ktrans*convterm(exp(theta),ttime[t]+t0,aifsettings),2);
     }
   return p;
 }

double space_space_update_theta2(double theta,double* kep, double ktrans, double vp, double t0, double* conc, double* time, double tau_epsilon, double* tau_theta, double* tau_theta2,double* tau_theta3, int x, int y, int z, int X, int Y, int Z, int T, double sigma, double* aifsettings, int spatial, int* img_mask)
{
//	if(setprint==1){printf("Start space_space_update_theta2() \n");}

  double theta_new=normal(theta,sigma);
  double logalpha = 0.0;
  if (spatial==3)
    {
      logalpha += space_log_fc_theta3(theta_new,kep,tau_epsilon,tau_theta,tau_theta2,tau_theta3,conc,time,t0,ktrans,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
      logalpha -= space_log_fc_theta3(theta,kep,tau_epsilon,tau_theta,tau_theta2,tau_theta3,conc,time,t0,ktrans,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
    }
  if (spatial==1)
    {
      logalpha += space_log_fc_theta2(theta_new,kep,tau_epsilon,tau_theta,tau_theta2,conc,time,t0,ktrans,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
      logalpha -= space_log_fc_theta2(theta,kep,tau_epsilon,tau_theta,tau_theta2,conc,time,t0,ktrans,vp,x,y,z,X,Y,Z,T,aifsettings, img_mask);
    }
  if (spatial==0)
    {
      logalpha += space_log_fc_theta(theta_new,tau_epsilon,tau_theta[ix(x,y,z,X,Y,Z)],conc,time,t0,ktrans,vp,x,y,z,X,Y,Z,T,aifsettings);
      logalpha -= space_log_fc_theta(theta,tau_epsilon,tau_theta[ix(x,y,z,X,Y,Z)],conc,time,t0,ktrans,vp,x,y,z,X,Y,Z,T,aifsettings);
    }

  if (exp(logalpha)>nulleins())
    {
     return theta_new;
    }
  else
    {
     return theta;
      
    }
}

void update_tau_global(double* tau_theta, double* tau_theta2,  double* tau_theta3, double a_theta, double a_theta3, double b_theta, double b_theta3, int* img_mask, double* kep, int X, int Y, int Z, int N, int spatial)
{
  // Update global precision, set same value for all local precisions
  double aa=a_theta;
  double bb=b_theta;
  int i, j;
  for (int x=1; x<X; x++)
    {
      for (int y=1; y<=Y; y++)
	{
	  for (int z=1; z<=Z; z++)
	    {
	      i=ix(x,y,z,X,Y,Z);
	      j=ix(x+1,y,z,X,Y,Z);
	      if (img_mask[i]==1)
		{
		  if (img_mask[j]==1)
		    {
		      aa+=0.5;
		      bb+=0.5*pow(log(kep[i])-log(kep[j]),2);
		    }
		}
	    }
	}
    }
  for (int x=1; x<=X; x++)
    {
      for (int y=1; y<Y; y++)
	{
	  for (int z=1; z<=Z; z++)
	    {
	      i=ix(x,y,z,X,Y,Z);
	      j=ix(x,y+1,z,X,Y,Z);
	      if (img_mask[i]==1)
		{
		  if (img_mask[j]==1)
		    {
		      aa+=0.5;
		      bb+=0.5*pow(log(kep[i])-log(kep[j]),2);
		    }
		}
	    }
	}
    }
  double tau=RNDGAM(aa,bb);
  //  Rprintf("\naa %f",aa);
  //  Rprintf(" bb %f\n",bb);

  for (int x=0; x<=N; x++)
    {
      tau_theta[x]=tau;
      tau_theta2[x]=tau;
    }

  // separate update for Z dimension
  if (spatial==3)
    {
      aa=a_theta;
      bb=b_theta;
      
      for (int x=1; x<=X; x++)
	{
	  for (int y=1; y<=Y; y++)
	    {
	      for (int z=1; z<Z; z++)
		{
		  i=ix(x,y,z,X,Y,Z);
		  j=ix(x,y,z+1,X,Y,Z);
		  if (img_mask[i]==1)
		    {
		      if (img_mask[j]==1)
			{
			  aa+=0.5;
			  bb+=0.5*pow(log(kep[i])-log(kep[j]),2);
			}
		    }
		}
	    }
	}
      tau=RNDGAM(aa,bb);

      for (int x=0; x<=N; x++)
	{
	  tau_theta3[x]=tau;
	}
    }
}

int space_tune(double* what,double acc, int i)
{
//	if(setprint==1){printf("Start space_tune() \n");}

  int count=0;

       if (acc<.3){what[i]=what[i]*.9;}
       if (acc<.2){what[i]=what[i]*.9;count=1;}
       if (acc<.1){what[i]=what[i]*.9;}
       if (acc<.02){what[i]=what[i]*.1;}
       if (acc>.5){what[i]=what[i]*1.1;}
       if (acc>.6){what[i]=what[i]*1.1;}
       if (acc>.7){what[i]=what[i]*1.1;count=1;}
       if (acc>.8){what[i]=what[i]*1.1;}
       if (acc>.9){what[i]=what[i]*1.1;}
       if (acc>.99){what[i]=what[i]*15;}

      return(count);
    
}

double update_tau_normal2(double tau, double a, double b, int* img_mask, double* theta, int X, int Y, int Z)
{
//	if(setprint==1){printf("Start update_tau_normal2() \n");}

  for (int x=1;x<=X; x++)
    {
      for (int y=1;y<=Y; y++)
	{
	  for (int z=1;z<=Z; z++)
	    {
	      if (img_mask[ix(x,y,z,X,Y,Z)])
		{
		  a++;
		  b+=0.5*pow(theta[ix(x,y,z,X,Y,Z)],2.0);
		}
	    }
	}
    }

  double temp = RNDGAM(a,b);
  return temp;
}

void update_tau_XY(double* tau_theta, double* tau_theta2,  double* tau_theta3, double a, double b, double a3, double b3, int* img_mask, double* theta, int X, int Y, int Z, int spatial)
{
//	if(setprint==1){printf("Start update_tau_XY() \n");}

  double bb;
  double aa;
  for (int x=1;x<X; x++)
  {
    for (int y=1;y<=Y; y++)
	{
	  for (int z=1;z<=Z; z++)
	  {
	    if (img_mask[ix(x,y,z,X,Y,Z)])
		{
		  if (img_mask[ix(x+1,y,z,X,Y,Z)])
		    {
		      if (x==1||x==(X-1))
			{
			  aa=a+1;
			  bb=b+0.5*pow(log(theta[ix(x+1,y,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0);
			}
		      else
			{
			  aa=a+2;
			  bb=b+pow(log(theta[ix(x+1,y,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0);
			}
		      tau_theta[ix(x,y,z,X,Y,Z)]=RNDGAM(aa,bb);
		    }
		}
	  }
  }
}

  for (int x=1;x<=X; x++)
    {
      for (int y=1;y<Y; y++)
	{
	  for (int z=1;z<=Z; z++)
	    {
	      if (img_mask[ix(x,y,z,X,Y,Z)])
		{
		  if (img_mask[ix(x,y+1,z,X,Y,Z)])
		    {
		      if (y==1||y==(Y-1))
			{
			  aa=a+1;		  
			  bb=b+0.5*pow(log(theta[ix(x,y+1,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0);
			}
		      else
			{
			  aa=a+2;		  
			  bb=b+pow(log(theta[ix(x,y+1,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0);
			}
		      tau_theta2[ix(x,y,z,X,Y,Z)]=RNDGAM(aa,bb);
		    }
		}
	    }
	}
    }

  if (spatial==3)
    {
      for (int x=1;x<=X; x++)
	{
	  for (int y=1;y<Y; y++)
	    {
	      for (int z=1;z<=Z; z++)
		{
		  if (img_mask[ix(x,y,z,X,Y,Z)])
		    {
		      if (img_mask[ix(x,y,z+1,X,Y,Z)])
			{
			  if (z==1||z==(Z-1))
			    {
			      aa=a3+1;		  
			      bb=b+0.5*pow(log(theta[ix(x,y,z+1,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0);
			    }
			  else
			    {
			      aa=a3+2;		  
			      bb=b+pow(log(theta[ix(x,y,z+1,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0);
			    }
			  tau_theta3[ix(x,y,z,X,Y,Z)]=RNDGAM(aa,bb);
			}
		    }
		}
	    }
	}
    }
  return;
}


// Update adaptive tau_theta for (tau_theta=tau_gamma), non-correct version, only 2D
void update_tau_XY2(double* tau_theta, double* tau_theta2,  double a, double b, double a3, double b3, int* img_mask, double* theta, double* tau_gamma, double* tau_gamma2, double* gamma, int X, int Y, int Z, int T)
{
//	if(setprint==1){printf("Start update_tau_XY2() \n");}

  double bb;
  double aa;
  for (int x=1;x<X; x++)
    {
      for (int y=1;y<=Y; y++)
	{
	  for (int z=1;z<=Z; z++)
	    {
	      if (img_mask[ix(x,y,z,X,Y,Z)])
		{
		  if (x==1||x==(X-1))
		    {
		      aa=a+2;
		      bb=b+0.5*pow(log(theta[ix(x+1,y,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0)+0.5*pow(log(gamma[ix(x+1,y,z,X,Y,Z)])-log(gamma[ix(x,y,z,X,Y,Z)]),2.0);
		    }
		  else
		    {
		      aa=a+4;
		      bb=b+pow(log(theta[ix(x+1,y,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0)+pow(log(gamma[ix(x+1,y,z,X,Y,Z)])-log(gamma[ix(x,y,z,X,Y,Z)]),2.0);
		    }
		  double temp=RNDGAM(aa,bb);
		  tau_gamma[ix(x,y,z,X,Y,Z)]=tau_theta[ix(x,y,z,X,Y,Z)]=temp;
		}
	    }
	}
    }      

  for (int x=1;x<=X; x++)
    {
      for (int y=1;y<Y; y++)
	{
	  for (int z=1;z<=Z; z++)
	    {
	      if (img_mask[ix(x,y,z,X,Y,Z)])
		{
		  if (y==1||y==(Y-1))
		    {
		      aa=a+2;		  
		      bb=b+0.5*pow(log(theta[ix(x,y+1,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0)+0.5*pow(log(gamma[ix(x,y+1,z,X,Y,Z)])-log(gamma[ix(x,y,z,X,Y,Z)]),2.0);
		    }
		  else
		    {
		      aa=a+4;		  
		      bb=b+pow(log(theta[ix(x,y+1,z,X,Y,Z)])-log(theta[ix(x,y,z,X,Y,Z)]),2.0)+pow(log(gamma[ix(x,y+1,z,X,Y,Z)])-log(gamma[ix(x,y,z,X,Y,Z)]),2.0);
		    }
		  double temp=RNDGAM(aa,bb);
		  tau_theta2[ix(x,y,z,X,Y,Z)]=tau_gamma2[ix(x,y,z,X,Y,Z)]=temp;
		}
	    }
	}
    }
  return;
}


double update_tau_epsilon(double a, double b, double* conc, int x, int y, int z, double vp, double ktrans, double kep, double t0, double* time, int X, int Y, int Z, int T, double* aifsettings)
{
//	if(setprint==1){printf("Start update_tau_epsilon() \n");}

  double bb=b;
  for (int t=1;t<=T; t++)
    {
      bb += .5*pow(conc[indx(x,y,z,t,X,Y,Z,T)]-extraterm(vp, time[t]+t0,aifsettings)-ktrans*convterm(kep, time[t]+t0,aifsettings),2.0);
    }

 /*
 if(setprint==1){
	  printf("aa= %f \n", a+0.5*(double)T);
	  printf("b= %f \n", b);
	  printf("bb= %f \n", bb);
	  printf("conc[1,1,2,22]: %f \n", conc[indx(1,1,2,22,X,Y,Z,T)]);
  }
*/

  return(RNDGAM(a+0.5*(double)T,bb));
}


double space_update_tau_epsilon1(double* tau, double aa, double bb, double* conc, int* img_mask, double* vp, double* ktrans, double* kep, double* t0, double* time, int X, int Y, int Z, int T, double* aifsettings)
{
 // if(setprint==1){printf("Start space_update_tau_epsilon1() \n");}
  
  double a=aa;
  double b=bb;

  //printf("pace_update_tau_epsilon1: apriori a %f, b %f \n", a,b);

  for (int x=1;x<=X; x++)
  {
      for (int y=1;y<=Y; y++)
	{
	  for (int z=1;z<=Z; z++)
	  {
	      if (img_mask[ix(x,y,z,X,Y,Z)])
	      {
	    	//printf("x,y,z: %d, %d, %d with ix(x,y,z,X,Y,Z): %d\n", x,y,z,ix(x,y,z,X,Y,Z));

	    	for (int t=1;t<=T; t++)
		    {
		      a += 0.5;
		      b += .5*pow(conc[indx(x,y,z,t,X,Y,Z,T)]-extraterm(vp[ix(x,y,z,X,Y,Z)], time[t]+t0[ix(x,y,z,X,Y,Z)],aifsettings)-ktrans[ix(x,y,z,X,Y,Z)]*convterm(kep[ix(x,y,z,X,Y,Z)], time[t]+t0[ix(x,y,z,X,Y,Z)],aifsettings),2.0);

		      //printf("(conc-fit)^2 %f \n",
		      //		  pow(conc[indx(x,y,z,t,X,Y,Z,T)]-extraterm(vp[ix(x,y,z,X,Y,Z)], time[t]+t0[ix(x,y,z,X,Y,Z)],aifsettings)-ktrans[ix(x,y,z,X,Y,Z)]*convterm(kep[ix(x,y,z,X,Y,Z)], time[t]+t0[ix(x,y,z,X,Y,Z)],aifsettings),2.0));


		    }
	      }
	   }
	  }
    }

  //printf("space_update_tau_epsilon1: aposteriori a %f, b %f \n", a,b);
  
  double tauneu=RNDGAM(a,b);
  if (!(tauneu>0)){tauneu=tau[1];}
  int N=X*Y*Z;
  for (int i=0; i<N; i++)
    {
      tau[i]=tauneu;
    }
  return(tauneu);

}

void dce_space(int* NRI,
	       double* conc,
	       int* img_mask,
	       int* dim,
	       double* ab_gamma, //5
	       double* ab_gamma3, 
	       double* ab_theta,
	       double* ab_theta3,
	       double* ab_vp,
	       double* ab_epsilon, //10
	       double* aif_settings,
	       int* settings,
	       double* time,
	       double* t0, //dim N!
	       double* tau_gamma, // dim N, 15
	       double* tau_theta, // dim N
	       double* ktrans, // dim N
	       double* kep, // dim N
	       double* vp, // dim N
	       double* tau_epsilon, // dim N, 20
	       double* ve, // dim N
	       double* sigmatheta, // dim N
	       double* sigmagamma, // dim N
	       double* sigmaeta,  // dim N
	       int* acc_gamma, // dim N, 25
	       int* acc_theta, // dim N
	       int* acc_eta, // dim N
	       int* acc,  // dim N
	       double* tau_gamma2,  // dim N
	       double* tau_theta2,  // dim N #30
	       double* tau_gamma3,  // dim N 
	       double* tau_theta3,  // dim N 
	       double* tau_eta2,  // dim N
	       double* tau_eta3,  // dim N
	       double* tau_eta,   // dim N #35
	       double* tau_vp,    // dim N
	       double* devi    // dim N
	       )  
{ 
	 GetRNGstate();

// dimensions
	 int X=dim[0];
	 int Y=dim[1];
	 int Z=dim[2];
	 int T=dim[3];
	 int N = X*Y*Z;
	 int N1 = N;

	 // parameters for priors
	 double a_theta=ab_theta[0];
	 double b_theta=ab_theta[1];
	 double a_gamma=ab_gamma[0];
	 double b_gamma=ab_gamma[1];
	 double a_vp=ab_vp[0];
	 double b_vp=ab_vp[1];
	 double a_theta3=ab_theta3[0];
	 double b_theta3=ab_theta3[1];
	 double a_gamma3=ab_gamma3[0];
	 double b_gamma3=ab_gamma3[1];
	 double a_epsilon=ab_epsilon[0];
	 double b_epsilon=ab_epsilon[1];

	 // setting parameters
	 int vpupdate=settings[0];
	 int spatial=settings[1];
	 int slice=settings[2];
	 int gemupdate=settings[3];
	 int uptauep=settings[4];
	 int respace_tunecycles=settings[5];
	 int space_tunepct=settings[6];

	 int nriters=NRI[0];
	 // int burnin=NRI[1] nicht nötig
	 int tuning=NRI[2];

	 double zaehler=0;

//	 printf("aifsettings model %f, a1 %f, m1 %f, a2 %f, m2 %f \n", aif_settings[0], aif_settings[1], aif_settings[2], aif_settings[3], aif_settings[4]);

	 double help;
	 int x,y,z;
	 
	 N1=0;
	 for (int x=1 ; x<=X  ; x++)
	 {
		 for (int y=1  ; y<=Y  ; y++)
		 {
			 for (int z=1; z<=Z; z++)
			 {
			   if(img_mask[ix(x,y,z,X,Y,Z)]==1){N1++;}
			 }
		  }
	   }
	   
	 //Rprint("%i pixels to analyse.\n",N1);
	 
	 int newvalue=0;
   
	 /************************************/
	 /* MCMC iterations                  */
	 /************************************/
	 //Rprintf("Starting iterations.\n",0);
	 double temp;
	 int respace_tune=0;
	 int sign=0;
	 for (int iter=1;iter<=nriters;iter++)
	 {
	   if (tuning>0)
     {
       if (iter==1)
       {
         Rprintf("Tuning:\nIteration 1");
       }
       if (fmod(iter,10.0)==0)
	     {
	   	    if(iter>100000){Rprintf("\b");}
     	    if(iter>10000){Rprintf("\b");}
     	    if(iter>1000){Rprintf("\b");}
	    	  if(iter>100){Rprintf("\b");}
	    	  if(iter>10){Rprintf("\b");}
	    	  Rprintf("\b%i",iter);
	    	}
     }

		 // iterate over all voxels
		 for (int riter=0;riter<2*N;riter++)
		 {
			 // randomly choose one voxel
			 int x=einsk(X);
			 int y=einsk(Y);
			 int z;
			 if (slice!=0)
			   {
			     z=slice;
    	   }
			 else
			   {
  		     z=einsk(Z);
  		   }
  
         if (img_mask[ix(x,y,z,X,Y,Z)]==1)
			   {
	         acc[ix(x,y,z,X,Y,Z)]=acc[ix(x,y,z,X,Y,Z)]+1;
			     
			     // update Ktrans
			     temp=space_update_gamma2(log(ktrans[ix(x,y,z,X,Y,Z)]), ktrans, kep[ix(x,y,z,X,Y,Z)], vp[ix(x,y,z,X,Y,Z)], tau_gamma, tau_gamma2,  tau_gamma3, tau_epsilon[ix(x,y,z,X,Y,Z)], conc, time, t0[ix(x,y,z,X,Y,Z)], x, y, z, X, Y, Z, T, sigmagamma[ix(x,y,z,X,Y,Z)], aif_settings, spatial, img_mask);
			     if (temp!=log(ktrans[ix(x,y,z,X,Y,Z)]))
			       {
				 acc_gamma[ix(x,y,z,X,Y,Z)]=acc_gamma[ix(x,y,z,X,Y,Z)]+1;
				 ktrans[ix(x,y,z,X,Y,Z)]=exp(temp);
			       }
			       
			     
			     // update kep
			     temp=space_space_update_theta2(log(kep[ix(x,y,z,X,Y,Z)]), kep, ktrans[ix(x,y,z,X,Y,Z)], vp[ix(x,y,z,X,Y,Z)], t0[ix(x,y,z,X,Y,Z)], conc, time, tau_epsilon[ix(x,y,z,X,Y,Z)], tau_theta, tau_theta2,tau_theta3, x, y, z, X, Y, Z, T, sigmatheta[ix(x,y,z,X,Y,Z)], aif_settings, spatial, img_mask);
			     if (temp!=log(kep[ix(x,y,z,X,Y,Z)]))
			       {
				 acc_theta[ix(x,y,z,X,Y,Z)]=acc_theta[ix(x,y,z,X,Y,Z)]+1;
				 kep[ix(x,y,z,X,Y,Z)]=exp(temp);
			       }
			     
			     
			     // update ve
			     ve[ix(x,y,z,X,Y,Z)]=ktrans[ix(x,y,z,X,Y,Z)]/kep[ix(x,y,z,X,Y,Z)];
			     
			     // update vp
			     if (vpupdate==1)
			       {
				 double v=T*tau_epsilon[ix(x,y,z,X,Y,Z)]+tau_vp[0];
				 double m=0.0;
				 for (int t=1; t<=T; t++)
				   {
				     m += tau_epsilon[ix(x,y,z,X,Y,Z)]*(conc[indx(x,y,z,t,X,Y,Z,T)]-ktrans[ix(x,y,z,X,Y,Z)]*convterm(kep[ix(x,y,z,X,Y,Z)],time[t]+t0[ix(x,y,z,X,Y,Z)],aif_settings));
				   }
				 m=m/v;
				 vp[ix(x,y,z,X,Y,Z)]=normal(m,v);
			       }
			     
			     if (vpupdate==2)
			       {
				 temp=space_update_eta2(logit(vp[ix(x,y,z,X,Y,Z)]), kep[ix(x,y,z,X,Y,Z)],  ktrans[ix(x,y,z,X,Y,Z)],  tau_eta,  tau_eta2,  tau_epsilon[ix(x,y,z,X,Y,Z)], conc, time,  t0[ix(x,y,z,X,Y,Z)], x,  y, z, X, Y, Z, T, sigmaeta[ix(x,y,z,X,Y,Z)], aif_settings);
				 if (temp!=logit(vp[ix(x,y,z,X,Y,Z)]))
				   {
				     acc_eta[ix(x,y,z,X,Y,Z)]=acc_eta[ix(x,y,z,X,Y,Z)]+1;
				     vp[ix(x,y,z,X,Y,Z)]=exp(temp)/exp(1+temp);
				   }
			       }
			     if (vpupdate==3)
			       {
				 temp=space_update_eta3(vp[ix(x,y,z,X,Y,Z)], kep[ix(x,y,z,X,Y,Z)],  ktrans[ix(x,y,z,X,Y,Z)],  tau_eta,  tau_eta2,  tau_epsilon[ix(x,y,z,X,Y,Z)], conc, time,  t0[ix(x,y,z,X,Y,Z)], x,  y, z, X, Y, Z, T, sigmaeta[ix(x,y,z,X,Y,Z)], aif_settings, a_vp, b_vp);
				 if (temp!=vp[ix(x,y,z,X,Y,Z)])
				   {
				     acc_eta[ix(x,y,z,X,Y,Z)]=acc_eta[ix(x,y,z,X,Y,Z)]+1;
				     vp[ix(x,y,z,X,Y,Z)]=temp;
				   }
			       }

			     if(spatial==0)
			     {
			    	 double taunew;
			    	 taunew=update_tau_epsilon(a_epsilon, b_epsilon, conc, x, y, z, vp[ix(x,y,z,X,Y,Z)], ktrans[ix(x,y,z,X,Y,Z)], kep[ix(x,y,z,X,Y,Z)], t0[ix(x,y,z,X,Y,Z)], time, X, Y, Z, T, aif_settings);
			    	 if(taunew >0)
			    	 {
			    		 tau_epsilon[ix(x,y,z,X,Y,Z)]=taunew;
			    	 }
			     }
			   } //endif img_mask[]==1
		 
		 } // end iterate over all voxels
					
		 // update tau_epsilon (inverse noise variance)
		 if (uptauep>=1)
		   {
		     space_update_tau_epsilon1(tau_epsilon, a_epsilon, b_epsilon, conc, img_mask, vp,ktrans,kep,t0,time, X, Y, Z, T, aif_settings);
		   }
		 
		 if (uptauep==3&&(spatial==1||spatial==3))
		   {
		     if (gemupdate==3)
		       {
			 // Update adaptive tau_theta for (tau_theta=tau_gamma), non-correct version, only 2D
			 update_tau_XY2(tau_gamma, tau_gamma2,  a_gamma, b_gamma, a_gamma3, b_gamma3, img_mask, ktrans,tau_theta, tau_theta2, kep, X, Y, Z, T);
		       }
		     else
		       {
			 // Update adaptive tau, non-correct version
			 update_tau_XY(tau_theta, tau_theta2, tau_theta3, a_theta, b_theta, a_theta3, b_theta3, img_mask, kep, X, Y, Z, T);
			 update_tau_XY(tau_gamma, tau_gamma2, tau_gamma3, a_gamma, b_gamma, a_gamma3, b_gamma3, img_mask, ktrans, X, Y, Z, T);
		       }
		   }
			
		 if (uptauep==2)
		   {
		     update_tau_global(tau_theta, tau_theta2, tau_theta3, a_theta, a_theta3, b_theta, b_theta3, img_mask, kep, X, Y, Z, N, spatial);
		     update_tau_global(tau_gamma, tau_gamma2, tau_gamma3, a_gamma, a_gamma3, b_gamma, b_gamma3, img_mask, kep, X, Y, Z, N, spatial);
       
       }

		 if (vpupdate==1&&spatial==0)
		   {
		     tau_vp[0]=update_tau_normal2(tau_vp[0], a_vp, b_vp,  img_mask, vp, X, Y, Z);
		   }

		 // tuning and acceptance rates
		 if (iter==tuning && respace_tune!=respace_tunecycles )
			 {
			   int count1 = 0;
			   int count2 = 0;
			   int count3 = 0;
			   int count4 = 0;
			   int count5 = 0;
			   int count6 = 0;
			   int tu=0;
			   for (int i=0; i<N; i++)
			     {
			       if (img_mask[i]==1)
				 {
				   count2 += acc_theta[i];
				   count3 += acc_gamma[i];
				   count4 += acc[i];
				   count6 += acc_eta[i];
				   if (gemupdate==0)
				     {
				       tu = space_tune(sigmatheta,(double)acc_theta[i]/(double)acc[i],i);
				       tu+= space_tune(sigmagamma,(double)acc_gamma[i]/(double)acc[i],i);
				       
				       if (vpupdate>=2)
					 {
					   tu+= space_tune(sigmaeta,(double)acc_eta[i]/(double)acc[i],i);
					 }
				     }
				   else
				     {
				       if (nulleins()>0.5)
					 {
					   tu = space_tune(sigmatheta,(double)acc_theta[i]/(double)acc[i],i);
					 }
				       else
					 {
					   tu =  space_tune(sigmagamma,(double)acc_gamma[i]/(double)acc[i],i);
					 }
				     }
				   if (tu>0)
				     {
				       count5++;
				     }
	
				 } //end if (img_mask[i])
			     } // end for


			   if (100*count5<N1*space_tunepct)
			     {
			       respace_tune++;
			       if(!(respace_tune==respace_tunecycles-1))
				 {
   			   Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bRe-tuning:\n");
				   iter=0;
				 }
			     }
			   else
			     {
     		   Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%i of ",count5);
           Rprintf("%i done.\n", N1);
			       iter=0;
			     }
			   for (int i=0; i<N; i++)
			     {
			       acc_theta[i]=0;
			       acc_eta[i]=0;
			       acc_gamma[i]=0;
			       acc[i]=0;			       
			     }
			 } // end if iter==tuning...


	 } // end for MCMC iterations
	 

		   PutRNGstate();

		   for (int x=1 ; x<=X  ; x++)
		   {
		     for (int y=1  ; y<=Y  ; y++)
			 {
			   for (int z=1; z<=Z; z++)
			   {
			       if(img_mask[ix(x,y,z,X,Y,Z)]==1)
				   {
			    	   double p = - T * log(tau_epsilon[ix(x,y,z,X,Y,Z)]);
			    	   for (int t=1; t<=T; t++)
			    	   {
			    		   p += tau_epsilon[ix(x,y,z,X,Y,Z)] * pow((conc[indx(x,y,z,t,X,Y,Z,T)]) - extraterm(vp[ix(x,y,z,X,Y,Z)], time[t], aif_settings) - ktrans[ix(x,y,z,X,Y,Z)] * convterm(kep[ix(x,y,z,X,Y,Z)], time[t], aif_settings), 2.0);
			    	   }
			    	   devi[ix(x,y,z,X,Y,Z)] = p;
				   }
			    }
			 }
		   }

		   return;
		   
} // end dce_space
