/*
##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
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
##
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "zufall.h"

double tune(double what,int acc,int n)
{
      if (acc<.3*n){what=what*.75;}
      if (acc<.2*n){what=what*.1;}
      if (acc<.1*n){what=what*.1;}
      if (acc<.02*n){what=what*.1;}
      if (acc>.6*n){what=what*1.1;}
      if (acc>.7*n){what=what*2;}
      if (acc>.8*n){what=what*5;}
      if (acc>.9*n){what=what*5;}
      if (acc>.99*n){what=what*15;}
         
      return(what);
    
}


double zahl(double mu, double sigma)
{
	static int gespeichert=0;
	static double alteZahl;

	if(gespeichert==0)
	{
		gespeichert=1;
		double X, Y, r2;

		do
		{
			X=reins();
			Y=reins();
			r2=X*X+Y*Y;		
		}while(r2>1.0);
	
		double norm=sqrt((-2.0*log(r2))/r2);
		alteZahl=Y*norm*sigma+mu;
		return X*norm*sigma+mu;
	}
	else
	{
		gespeichert=0;
		return alteZahl;
	}
}


double aif(double t, double* settings)
{
  double aifmodel=settings[0];
  double a1=settings[1];
  double m1=settings[2];
  double a2=settings[3];
  double m2=settings[4];
  double aif = 0.0;
  if (aifmodel==0)
    {
      aif = (a1*exp(-m1*t)+a2*exp(-m2*t));
    }
  if (aifmodel==1)
    {
      aif=a1*t*exp(-m1*t)+a2*(exp(-m2*t)-exp(-m1*t));
    }
  return(aif);
}

double extraterm(double vp,double time,double* settings)
{
      return(vp*aif(time,settings));
}


double convterm(double kep, double t, double* settings)
{
  double aifmodel=settings[0];
  double a1=settings[1];
  double m1=settings[2];
  double a2=settings[3];
  double m2=settings[4];
  double bla=0.0;

  if (aifmodel==0)
    {
      bla = a1*(exp(-m1*t)-exp(-kep*t))/(kep-m1)+a2*(exp(-m2*t)-exp(-kep*t))/(kep-m2);
    }
  if (aifmodel==1)
    {
      bla = a1*kep* (t*exp(-m1*t) - (exp(-m1*t)-exp(-kep*t))/(kep-m1))/(kep-m1);
      bla += a2*kep*( (exp(-m2*t)-exp(-kep*t))/(kep-m2) - (exp(-m1*t)-exp(-kep*t))/(kep-m1) );
    }
  return bla;
}

 
double log_fc_gamma(double gamma, double tau_epsilon, double tau_gamma, double* conc, double* time, double kep, double vp, int T, double* settings)
 {
   int t;
   double p = 0.0;
   p -= 0.5*tau_gamma*gamma*gamma;
   for (t=0; t<T; t++)
     {
       p -= 0.5*tau_epsilon*pow((conc[t])-extraterm(vp,time[t],settings)-exp(gamma)*convterm(kep,time[t],settings),2);
     }
   return p;
 }

double update_gamma2(double gamma, double kep, double vp, double tau_gamma, double tau_epsilon, double* conc, double* time, double sigma, int T, double* aif_settings){

  double gamma_new = normal(gamma,sigma);
  double logalpha = 0.0;
  logalpha +=  log_fc_gamma(gamma_new,tau_epsilon,tau_gamma,conc,time,kep,vp,T,aif_settings);
  logalpha -= log_fc_gamma(gamma,tau_epsilon,tau_gamma,conc,time,kep,vp,T,aif_settings);


  if (exp(logalpha)>nulleins())
    {
      return gamma_new;
    }
  else
    {
      return gamma;
    }
}

double log_fc_theta(double theta, double tau_epsilon, double tau_theta, double* conc, double* time,double ktrans, double vp, int T, double* settings)
 {
   double p = 0.0;
   int t;
   p -= 0.5*tau_theta*(theta+1)*(theta+1);
   for (t=0; t<T; t++)
     {
       p -= 0.5*tau_epsilon*pow((conc[t])-extraterm(vp,time[t],settings)-ktrans*convterm(ktrans/exp(theta),time[t],settings),2);
     }
   return p;
 }
double update_theta2(double theta, double ktrans, double vp, double* conc, double* time, double tau_epsilon, double tau_theta, double sigma, int T, double* aif_settings)
{
  double theta_new=normal(theta,sigma);
  double logalpha = 0.0;
  logalpha += log_fc_theta(theta_new,tau_epsilon,tau_theta,conc,time,ktrans,vp,T,aif_settings);
  logalpha -= log_fc_theta(theta,tau_epsilon,tau_theta,conc,time,ktrans,vp,T,aif_settings);
 
  if (exp(logalpha)>nulleins())
    {
     return theta_new;
    }
  else
    {
     return theta;
    }
}
double log_fc_eta3(double eta, double tau_epsilon, double a_vp, double b_vp, double* conc, double* time, double kep, double ktrans, int T, double* settings)
 {
   double p = 0.0;
   double a,b;
   int t;
   for (t=0; t<T; t++)
     {
       a = aif(time[t],settings);
       b=conc[t]-ktrans*convterm(kep,time[t],settings);
       p -= 0.5*tau_epsilon*(a*eta-b)*(a*eta-b);
     }
   p += (a_vp-1)*log(eta)+(b_vp-1)*log(1-eta);
   return p;
 }
double update_eta3(double vp, double kep, double ktrans, double a_vp, double b_vp, double tau_epsilon, double* conc, double* time, double sigma, int T,double* aif_settings){
  double eta_new = -1;
  while (eta_new>1 || eta_new<0)
    {
      eta_new=normal(vp,sigma);
    }
  double logalpha = 0.0;
  logalpha +=  log_fc_eta3(eta_new,tau_epsilon,a_vp, b_vp,conc,time,kep,ktrans,T,aif_settings);
  logalpha -= log_fc_eta3(vp,tau_epsilon,a_vp, b_vp,conc,time,kep,ktrans,T,aif_settings);
  if (exp(logalpha)>nulleins())
    {
      return eta_new;
    }
  else
    {
      return vp;
    }
}
double update_tau_epsilon1(double tau, double aa, double bb, double* conc, double vp, double ktrans, double kep, double* time, int T, double* settings)
{
  int i;
  for (i=0;i<T; i++)
    {
      aa+=0.5;
      bb += .5*pow(conc[i]-extraterm(vp, time[i],settings)-ktrans*convterm(kep, time[i],settings),2.0);
    }

  tau=RNDGAM(aa,bb);
  return(tau);
}





void dce_bayes_run_single(int* NRI, 
	 double* conc,
	 double* tau_gamma, double* tau_theta,
	 double* ab_vp,
	 double* ab_epsilon, 
	 double* aif_settings, int* settings, double* time, int* T,
	 double *ktrans_trace, double* kep_trace, double* vp_trace, double* tau_epsilon_trace)
	 
{
  GetRNGstate();
  int iter,tu;
  double temp;
  double ktrans=.5;
  double kep=1;
  double vp=0.0;
  if (settings[0]==1)
    vp=ab_vp[0]/ab_vp[1];
  double tau_epsilon=ab_epsilon[0]/ab_epsilon[1];
  int sample=-1;
  double sigmagamma=1;
  double sigmatheta=1;
  double sigmaeta=1;
  int acc_gamma=0;
  int acc_theta=0;
  int acc_eta=0;
  iter=0;

  while (iter<NRI[0])
    { 
      iter++;
     temp=update_gamma2(log(ktrans), kep, vp,  tau_gamma[0], tau_epsilon, conc, time, sigmagamma, T[0], aif_settings);
     
    if (temp!=log(ktrans))
	{
	  acc_gamma++;
	  ktrans=exp(temp);
	}
      
      temp=update_theta2(log(ktrans/kep), ktrans, vp,  conc, time, tau_epsilon, tau_theta[0], sigmatheta, T[0], aif_settings);
      if (temp!=log(ktrans/kep))
	{
	  acc_theta++;
	  kep=ktrans/exp(temp);
	}
      
  
      if (settings[0]==1)
	{
	  temp=update_eta3(vp, kep, ktrans, ab_vp[0], ab_vp[1], tau_epsilon, conc, time, sigmaeta, T[0], aif_settings);
	  if (temp!=vp)
	    {
	      acc_eta++;
	      vp=temp;
	    }
	}
      
      tau_epsilon=update_tau_epsilon1(tau_epsilon, ab_epsilon[0], ab_epsilon[1], conc, vp,ktrans,kep,time,T[0], aif_settings);
    

      if (iter==NRI[3])
	{

	  tu=0;
	     temp=tune(sigmagamma,acc_gamma,NRI[3]);
	     if (sigmagamma!=temp)
	       {
		 sigmagamma=temp;
		 tu=1;
	       }
	  temp=tune(sigmatheta,acc_theta,NRI[3]);
	  if (sigmatheta!=temp)
	    {
	      sigmatheta=temp;
	      tu=1;
	    }
	  if(settings[0]==1)
	    {
	      temp=tune(sigmaeta,acc_eta,NRI[3]);
	      if (sigmaeta!=temp)
		{
		  sigmaeta=temp;
		  tu=1;
		}
	    }
	  if (tu!=0)
	    {
	      iter=0;
	      acc_gamma=0;
	      acc_theta=0;
	      acc_eta=0;
	    }
	  else
	    {
	      sample=0;
	    }
	}

      if (iter>NRI[2] && fmod(iter,NRI[1])==0 && sample>=0)
       {
	 ktrans_trace[sample]=ktrans;
	 kep_trace[sample]=kep;
	 vp_trace[sample]=vp;
	 tau_epsilon_trace[sample]=tau_epsilon;
	 sample++;
       }


    }
      
 
  PutRNGstate();

}
