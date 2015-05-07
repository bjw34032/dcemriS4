/*
##
## Copyright (c) 2012, Julia Sommer, Volker Schmid
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



/******************************/
/* Global Variables			  */
/******************************/
//int printerlevel = 0;

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "zufall.h"
#include "baymri.h"
#include "space.h"

/************************************************************************************************/
double log_likelihood(double tau_epsilon, double ktrans, double kep, double ktrans2, double kep2,
						double* conc, double* time, double vp,
						int x, int y, int z, int X, int Y, int Z, int T, double* aifsettings)
{
   double p = 0.0;
   p +=  T * log(tau_epsilon);
   for (int t=1; t<=T; t++)
   {
	   p -= tau_epsilon * pow((conc[indx(x,y,z,t,X,Y,Z,T)])
			   - extraterm(vp, time[t], aifsettings)
			   - ktrans * convterm(kep, time[t], aifsettings)
			   - ktrans2 * convterm(kep2, time[t], aifsettings), 2.0);
   }

   p = 0.5* p;
   return p;
}

/************************************************************************************************************/
double log_priori_gamma(double gamma, double* ktrans, double* tau_gamma, double* tau_gamma2,
						int x, int y, int z, int X, int Y, int Z, int* img_mask)
{
   double p = 0.0;
   if (x!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma[ix(x-1,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x-1,y,z,X,Y,Z)]),2.0);}}
   if (x!=X){if(img_mask[ix(x+1,y,z,X,Y,Z)]){p -= 0.5*tau_gamma[ix(x,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x+1,y,z,X,Y,Z)]),2.0);}}
   if (y!=1){if(img_mask[ix(x,y-1,z,X,Y,Z)]){p -= 0.5*tau_gamma2[ix(x,y-1,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y-1,z,X,Y,Z)]),2.0);}}
   if (y!=Y){if(img_mask[ix(x,y+1,z,X,Y,Z)]){p -= 0.5*tau_gamma2[ix(x,y,z,X,Y,Z)]*pow(gamma-log(ktrans[ix(x,y+1,z,X,Y,Z)]),2.0);}}
   return p;
}

/************************************************************************************************************/
double log_priori_theta(double theta, double* kep, double* tau_theta, double* tau_theta2,
						int x, int y, int z, int X, int Y, int Z, int* img_mask)
 {
   double p = 0.0;
   if (x!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_theta[ix(x-1,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x-1,y,z,X,Y,Z)]),2.0);}}
   if (x!=X){if(img_mask[ix(x+1,y,z,X,Y,Z)]){p -= 0.5*tau_theta[ix(x,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x+1,y,z,X,Y,Z)]),2.0);}}
   if (y!=1){if(img_mask[ix(x,y-1,z,X,Y,Z)]){p -= 0.5*tau_theta2[ix(x,y-1,z,X,Y,Z)]*pow(theta-log(kep[ix(x,y-1,z,X,Y,Z)]),2.0);}}
   if (y!=Y){if(img_mask[ix(x,y+1,z,X,Y,Z)]){p -= 0.5*tau_theta2[ix(x,y,z,X,Y,Z)]*pow(theta-log(kep[ix(x,y+1,z,X,Y,Z)]),2.0);}}

   return p;
 }


/************************************************************************************************************/
double log_priori_beta(double beta, double* kep, double* ktrans, double* tau_beta, double* tau_beta2,
						int x, int y, int z, int X, int Y, int Z, int* img_mask)
 {
   double p = 0.0;
   if (x!=1){if(img_mask[ix(x-1,y,z,X,Y,Z)]){p -= 0.5*tau_beta[ix(x-1,y,z,X,Y,Z)]*pow(beta-log(ktrans[ix(x-1,y,z,X,Y,Z)]/kep[ix(x-1,y,z,X,Y,Z)]),2.0);}}
   if (x!=X){if(img_mask[ix(x+1,y,z,X,Y,Z)]){p -= 0.5*tau_beta[ix(x,y,z,X,Y,Z)]*pow(beta-log(ktrans[ix(x+1,y,z,X,Y,Z)]/kep[ix(x+1,y,z,X,Y,Z)]),2.0);}}
   if (y!=1){if(img_mask[ix(x,y-1,z,X,Y,Z)]){p -= 0.5*tau_beta2[ix(x,y-1,z,X,Y,Z)]*pow(beta-log(ktrans[ix(x,y-1,z,X,Y,Z)]/kep[ix(x,y-1,z,X,Y,Z)]),2.0);}}
   if (y!=Y){if(img_mask[ix(x,y+1,z,X,Y,Z)]){p -= 0.5*tau_beta2[ix(x,y,z,X,Y,Z)]*pow(beta-log(ktrans[ix(x,y+1,z,X,Y,Z)]/kep[ix(x,y+1,z,X,Y,Z)]),2.0);}}

   return p;
 }


/************************************************************************************************************/
double log_priori_nonspatial(double theta, double tau_theta)
 {
	 double p = 0.0;
	 p -= 0.5*tau_theta*theta*theta;
	 return p;
 }


/************************************************************************************************************/
double update_gamma(double gamma, double* ktrans, double kep, double ktrans2, double kep2, double vp,
					double* tau_gamma, double* tau_gamma2, double* tau_gamma3,double tau_epsilon,
					double* conc, double* time, double t0, int x, int y, int z, int X, int Y, int Z, int T,
					double sigma, double* aifsettings, int spatial, int* img_mask)
{

//	if(printerlevel==1){printf("Start update_gamma() \n");}

	double gamma_new = normal(gamma,sigma);
	double logalpha = 0.0;

    if (spatial==1)
    {
    	logalpha += log_likelihood(tau_epsilon, exp(gamma_new), kep, ktrans2, kep2,
									conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

    	logalpha += log_priori_gamma(gamma_new,ktrans,tau_gamma,tau_gamma2,x,y,z,X,Y,Z,img_mask);

    	logalpha -= log_likelihood(tau_epsilon, exp(gamma), kep, ktrans2, kep2,
    										conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

    	logalpha -= log_priori_gamma(gamma,ktrans,tau_gamma,tau_gamma2,x,y,z,X,Y,Z,img_mask);
    }

    if (spatial==0)
    {
    	logalpha += log_likelihood(tau_epsilon, exp(gamma_new), kep, ktrans2, kep2,
									conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

    	// non-spatial prior
    	logalpha += log_priori_nonspatial(gamma_new, tau_gamma[ix(x,y,z,X,Y,Z)]);


    	logalpha -= log_likelihood(tau_epsilon, exp(gamma), kep, ktrans2, kep2,
    										conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

    	// non-spatial prior
    	logalpha -= log_priori_nonspatial(gamma, tau_gamma[ix(x,y,z,X,Y,Z)]);
    }

	if (exp(logalpha)>nulleins())
    {
		//printf("new ktrans accepted \n");
      return gamma_new;
    }
	else
    {
		//printf("keep old ktrans \n");
	  return gamma;
    }
}

/************************************************************************************************************/
double update_gamma3(double gamma, double* ktrans, double kep, double ktrans2, double kep2, double vp,
					double* tau_gamma, double* tau_gamma2, double* tau_gamma3,double tau_epsilon,
					double* conc, double* time, double t0, int x, int y, int z, int X, int Y, int Z, int T,
					double sigma, double* aifsettings, int spatial, int* img_mask)
{

//	if(printerlevel==1){printf("Start update_gamma3() \n");}

	double gamma_new = normal(gamma,sigma);
	double logalpha = 0.0;

    if (spatial==1)
    {
    	logalpha += log_likelihood(tau_epsilon, exp(gamma_new), kep, ktrans2, kep2,
									conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );
    	//printf("logalpha log_likeli(gamma_new) %f\n", logalpha);

    	logalpha += log_priori_gamma(gamma_new,ktrans,tau_gamma,tau_gamma2,x,y,z,X,Y,Z,img_mask);
    	//printf("logalpha log_prior(gamma_new) %f\n", logalpha);

    	logalpha -= log_likelihood(tau_epsilon, exp(gamma), kep, ktrans2, kep2,
    										conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );
    	//printf("logalpha log_likeli(gamma) %f\n", logalpha);

    	logalpha -= log_priori_gamma(gamma,ktrans,tau_gamma,tau_gamma2,x,y,z,X,Y,Z,img_mask);
    }

    if (spatial==0)
    {

    	double help= 0.0;
    	logalpha += log_likelihood(tau_epsilon, exp(gamma_new), kep, ktrans2, kep2,
									conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

    	// non-spatial prior
    	help = gamma_new-log(5); // Priori ktrans2 with expectation 5
    	logalpha += log_priori_nonspatial(help, tau_gamma[ix(x,y,z,X,Y,Z)]);


    	logalpha -= log_likelihood(tau_epsilon, exp(gamma), kep, ktrans2, kep2,
    										conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

    	// non-spatial prior
    	help = gamma-log(5);
    	logalpha -= log_priori_nonspatial(help, tau_gamma[ix(x,y,z,X,Y,Z)]);
    }

	if (exp(logalpha)>nulleins())
    {
		//printf("new ktrans accepted \n");
      return gamma_new;
    }
	else
    {
		//printf("keep old ktrans \n");
	  return gamma;
    }
}

/************************************************************************************************************/
double update_theta(double theta,double* kep, double ktrans, double kep2, double ktrans2, double vp,
		 double* tau_theta, double* tau_theta2,double* tau_theta3, double tau_epsilon,
		 double* conc, double* time, double t0,
		int x, int y, int z, int X, int Y, int Z, int T, double sigma, double* aifsettings, int spatial, int* img_mask)
{
//	if(printerlevel==1){printf("Start update_theta() \n");}

	double theta_new = normal(theta,sigma);
	double logalpha = 0.0;

	if (spatial==1)
	{
		logalpha += log_likelihood(tau_epsilon, ktrans, exp(theta_new), ktrans2, kep2,
										conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

		//printf("log_likeli(theta_new) %f\n", log_likelihood(tau_epsilon, ktrans, exp(theta_new), ktrans2, kep2,conc, time, vp, x,y,z,X,Y,Z,T, aifsettings ));

		logalpha += log_priori_theta(theta_new, kep, tau_theta, tau_theta2, x,y,z,X,Y,Z,img_mask);

		//printf("logalpha log_priori(theta_new) %f\n", logalpha);

		logalpha -= log_likelihood(tau_epsilon, ktrans, exp(theta), ktrans2, kep2,
										conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

		//printf("log_likeli(theta) %f\n", log_likelihood(tau_epsilon, ktrans, exp(theta), ktrans2, kep2,conc, time, vp, x,y,z,X,Y,Z,T, aifsettings ));

		logalpha -= log_priori_theta(theta, kep, tau_theta, tau_theta2, x,y,z,X,Y,Z,img_mask);
	}

	if (spatial==0)
	{
		logalpha += log_likelihood(tau_epsilon, ktrans, exp(theta_new), ktrans2, kep2,
											conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );
		// non-spatial prior
		logalpha += log_priori_nonspatial(theta_new, tau_theta[ix(x,y,z,X,Y,Z)]);


		logalpha -= log_likelihood(tau_epsilon, ktrans, exp(theta), ktrans2, kep2,
											conc, time, vp, x,y,z,X,Y,Z,T, aifsettings );

		// non-spatial prior
		logalpha -= log_priori_nonspatial(theta, tau_theta[ix(x,y,z,X,Y,Z)]);
	}

	if (exp(logalpha)>nulleins())
	{
		//printf("new kep accepted \n");
		return theta_new;
	}
	else
	{
		//printf("return old kep \n");
		return theta;
	}
}

/************************************************************************************************************/
double log_fc_vp(double vp, double tau_epsilon, double tau_vp, double* conc, double* time,
				double t0, double kep, double ktrans, double kep2, double ktrans2,
				int x, int y, int z, int X, int Y, int Z, int T,
				double* aifsettings, double a_vp, double b_vp)
{
   double p = 0.0;
   double a,b;
   for (int t=1; t<=T; t++)
   {
       a = aif(time[t]+t0,aifsettings);
       b = conc[indx(x,y,z,t,X,Y,Z,T)]
                -ktrans*convterm(kep,time[t]+t0,aifsettings)
                -ktrans2*convterm(kep2,time[t]+t0,aifsettings);
       p -= 0.5*tau_epsilon*(a*vp-b)*(a*vp-b);
   }
   p += (a_vp-1)*log(vp)+(b_vp-1)*log(1-vp); // log prior, beta prior on vp
   return p;
}

/************************************************************************************************************/
double update_vp(double vp, double kep, double ktrans, double kep2, double ktrans2,
		double* tau_vp, double* tau_vp2, double tau_epsilon,
		double* conc, double* time, double t0,
		int x, int y, int z, int X, int Y, int Z, int T, double sigma,
		double* aifsettings,double a_vp, double b_vp)
{
//	if(printerlevel==1){printf("Start update_vp() \n");}

	double vp_new = -1;
	while (vp_new>1 || vp_new<0)
    {
		vp_new=normal(vp,sigma);
    }
	double logalpha = 0.0;

	logalpha += log_fc_vp(vp_new,tau_epsilon,tau_vp[ix(x,y,z,X,Y,Z)],
				conc,time,t0,kep,ktrans,kep2,ktrans2,x,y,z,X,Y,Z,T,aifsettings,a_vp,b_vp);
	logalpha -= log_fc_vp(vp,tau_epsilon,tau_vp[ix(x,y,z,X,Y,Z)],
				conc,time,t0,kep,ktrans,kep2,ktrans2,x,y,z,X,Y,Z,T,aifsettings,a_vp,b_vp);

	if (exp(logalpha)>nulleins())
    {
      return vp_new;
    }
	else
    {
      return vp;
    }
}

/************************************************************************************************************/
double update_tau_epsilon_voxelwise(double a, double b, double* conc, int x, int y, int z, double vp,
						double ktrans, double kep, double ktrans2, double kep2, double t0, double* time,
						int X, int Y, int Z, int T, double* aifsettings)
{
//	  if(printerlevel==1){printf("Start update_tau_epsilon_voxelwise() \n");}

	  double bb=b;
	  for (int t=1;t<=T; t++)
	  {
		  bb += .5*pow(conc[indx(x,y,z,t,X,Y,Z,T)]
						-extraterm(vp, time[t]+t0,aifsettings)
						-ktrans*convterm(kep, time[t]+t0,aifsettings)
						-ktrans2*convterm(kep2, time[t]+t0,aifsettings),2.0);
	  }

/*
	  if(printerlevel==1){
		  printf("aa= %f \n", a+0.5*(double)T);
		  printf("b= %f \n", b);
		  printf("ktrans2= %f \n", ktrans2);
		  printf("kep2= %f \n", kep2);
		  printf("bb= %f \n", bb);
		  printf("conc[1,1,2,22]-ktrans2*convterm...: %f \n", conc[indx(1,1,2,22,X,Y,Z,T)]-ktrans2*convterm(kep2, time[21]+t0,aifsettings));
		  printf("ktrans2*convterm...: %f \n", ktrans2*convterm(kep2, time[21]+t0,aifsettings));

	  }
*/

	  return(RNDGAM(a+0.5*(double)T,bb));
}

/************************************************************************************************************/
double update_tau_epsilon_all(double* tau, double aa, double bb, double* conc, int* img_mask,
		double* vp, double* ktrans, double* kep, double* ktrans2, double* kep2,
		double* t0, double* time, int X, int Y, int Z, int T, double* aifsettings)
{
//	if(printerlevel==1){printf("Start update_tau_epsilon_all() \n");}

	double a=aa;
	double b=bb;

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
						b += .5*pow(conc[indx(x,y,z,t,X,Y,Z,T)]
						      - extraterm(vp[ix(x,y,z,X,Y,Z)], time[t]+t0[ix(x,y,z,X,Y,Z)],aifsettings)
						      - ktrans[ix(x,y,z,X,Y,Z)]*convterm(kep[ix(x,y,z,X,Y,Z)], time[t]+t0[ix(x,y,z,X,Y,Z)],aifsettings)
						      - ktrans2[ix(x,y,z,X,Y,Z)]*convterm(kep2[ix(x,y,z,X,Y,Z)], time[t]+t0[ix(x,y,z,X,Y,Z)],aifsettings),2.0);

					}
				}
			}
		}
    }

	//printf("update_tau_epsilon_all: aposteriori a %f, b %f \n", a,b);

	double tauneu=RNDGAM(a,b);
	if (!(tauneu>0))
	{
		tauneu=tau[1];
	}
	int N=X*Y*Z;
	for (int i=0; i<N; i++)
	{
	   tau[i]=tauneu;
	}
	return(tauneu);

}

/************************************************************************************************************/
void update_tau_theta(double* tau_theta, double* tau_theta2,  double* tau_theta3,
						double a_theta, double a_theta3, double b_theta, double b_theta3,
						int* img_mask, double* kep, int X, int Y, int Z, int N, int spatial)
{
	// Update global precision,  (same value for all local precisions)
	double aa=a_theta;
	double bb=b_theta;
	int i, j;

	// neighbors in x-direction
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
	// neighbors in y-direction
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

	// draw new tau and save it (same tau for x- and y-direction)
	double tau=RNDGAM(aa,bb);
	for (int x=0; x<=N; x++)
    {
		tau_theta[x]=tau;
		tau_theta2[x]=tau;
    }

	// neighbors in z-direction
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
		// draw new tau and save it (separate tau for z-direction)
     	tau=RNDGAM(aa,bb);

     	for (int x=0; x<=N; x++)
     	{
     		tau_theta3[x]=tau;
     	}
    }
}

/****************************************************************************************/
/* main                                                                                 */
/****************************************************************************************/
void dce_space_2comp(
	// input variables
		   int* NRI,
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
	       double* time, // dim T+1
	       double* t0, //dim N!
	// input-output variables
	       double* ktrans, // dim N //15
	       double* kep, // dim N
	       double* ktrans2, // dim N
	       double* kep2, // dim N
	       double* vp, // dim N
	       double* tau_epsilon, // dim N, //20
	       // precisions x-y-direction
	       double* tau_gamma, // dim N,
	       double* tau_theta, // dim N
	       double* tau_gamma_2comp,  // dim N
		   double* tau_theta_2comp,  // dim N
		   double* tau_vp,   // dim N //25
		   // precisions z-direction
		   double* tau_gamma3,  // dim N
		   double* tau_theta3,  // dim N
		   double* tau_gamma3_2comp,  // dim N
		   double* tau_theta3_2comp,  // dim N
		   double* tau_vp3,  // dim N //30
		   // proposal variances
	       double* sigmatheta, // dim N
	       double* sigmagamma, // dim N
	       double* sigmatheta_2comp, // dim N
	       double* sigmagamma_2comp, // dim N
	       double* sigmavp,  // dim N  //35
	 // output variables
	       // acceptance rates
	       int* acc_gamma, // dim N,
	       int* acc_theta, // dim N
	       int* acc_gamma_2comp, // dim N,
	       int* acc_theta_2comp, // dim N
	       int* acc_vp, // dim N //40
	       int* acc,  // dim N
	       // deviance
	       double* devi,    // dim N
	       // fit
	       double* fit    // dim N*T
	       )
{ 
	 GetRNGstate();
     //   Rprintf("%f - ",ab_gamma[0]);
    //    Rprintf("%f\n",ab_gamma[1]);

	// Rprintf("Start C-procedure dce_space_2comp\n",0);

	 // dimensions
	 int X=dim[0];
	 int Y=dim[1];
	 int Z=dim[2];
	 int T=dim[3];
	 int N = X*Y*Z+1;
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
//	 int regularize=settings[7];

	 // precision in x-y-direction is the same
	 double tau_gamma2[N];
	 double tau_gamma2_2comp[N];
	 double tau_theta2[N];
	 double tau_theta2_2comp[N];
	 double tau_vp2[N];
	 for(int i=0; i<N; i++)
	 {
		 tau_gamma2[i] = tau_gamma[i];
		 tau_gamma2_2comp[i] = tau_gamma_2comp[i];
		 tau_theta2[i] = tau_theta[i];
		 tau_theta2_2comp[i] = tau_theta_2comp[i];
		 tau_vp2[i] = tau_vp[i];
	 }

	 // volumes
	 double ve[N];
	 double ve2[N];
	 for(int i=0; i<N; i++)
	 {
		 if(img_mask[i]==1)
		 {
			 ve[i] = ktrans[i]/kep[i];
			 ve2[i] = ktrans2[i]/kep2[i];
		 }
	 }


	 int nriters=NRI[0];
	 int tuning=NRI[2];
	 double zaehler=0;
	 int x,y,z;
	 double help;

/*	 printf("Starting values parameter: \n");
	 printf("ktrans[1,1,2] %f \n", ktrans[ix(1,1,2,X,Y,Z)]);
	 printf("ktrans2[1,1,2] %f \n", ktrans2[ix(1,1,2,X,Y,Z)]);
	 printf("kep[1,1,2] %f \n", kep[ix(1,1,2,X,Y,Z)]);
	 printf("kep2[1,1,2] %f \n", kep2[ix(1,1,2,X,Y,Z)]);

	 printf("tau_epsilon[1,1,2] %f \n", tau_epsilon[ix(1,1,2,X,Y,Z)]);
	 printf("tau_gamma[1,1,2] %f \n", tau_gamma[ix(1,1,2,X,Y,Z)]);
	 printf("tau_theta[1,1,2] %f \n", tau_theta[ix(1,1,2,X,Y,Z)]);
	 printf("tau_gamma_2comp[1,1,2] %f \n", tau_gamma_2comp[ix(1,1,2,X,Y,Z)]);
	 printf("tau_theta_2comp[1,1,2] %f \n", tau_theta_2comp[ix(1,1,2,X,Y,Z)]);

	 printf("sigmagamma[1,1,2] %f \n", sigmagamma[ix(1,1,2,X,Y,Z)]);
	 printf("sigmagamma_2comp[1,1,2] %f \n", sigmagamma_2comp[ix(1,1,2,X,Y,Z)]);
	 printf("sigmatheta[1,1,2] %f \n", sigmatheta[ix(1,1,2,X,Y,Z)]);
	 printf("sigmatheta_2comp[1,1,2] %f \n", sigmatheta_2comp[ix(1,1,2,X,Y,Z)]);*/

//	 printf("aifsettings model %f, a1 %f, m1 %f, a2 %f, m2 %f \n", aif_settings[0], aif_settings[1], aif_settings[2], aif_settings[3], aif_settings[4]);
//	 printf("settings: spatial %d, uptauep %d  \n", spatial, uptauep);

	 
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

	 Rprintf("Starting iterations.\n",0);
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
	   		   if(iter>10){Rprintf("\b");}
     		   if(iter>100){Rprintf("\b");}
     		   if(iter>1000){Rprintf("\b");}
       	   if(iter>10000){Rprintf("\b");}
       	   if(iter>100000){Rprintf("\b");}
	    	   Rprintf("\b%i",iter);
	   	   }
	    	}

       // Rprintf("%f - ",tau_theta[0]);
//        Rprintf("%f\n",tau_gamma[0]);

    // iterate over all voxels
		 for (int riter=0;riter<2*N;riter++)
		 {
			 // randomly choose one voxel
			 int x=einsk(X);
			 int y=einsk(Y);
			 int z;
			 if (slice!=0) // slice fixed
			 {
				 z=slice;
			 }
			 else // random slice
			 {
			     z=einsk(Z);
			 }
			 if (img_mask[ix(x,y,z,X,Y,Z)]==1)
			 {
			     acc[ix(x,y,z,X,Y,Z)]=acc[ix(x,y,z,X,Y,Z)]+1;
			     
			     // update kep
			     temp = update_theta(log(kep[ix(x,y,z,X,Y,Z)]), kep, ktrans[ix(x,y,z,X,Y,Z)],kep2[ix(x,y,z,X,Y,Z)],ktrans2[ix(x,y,z,X,Y,Z)],
							 vp[ix(x,y,z,X,Y,Z)], tau_theta, tau_theta2,tau_theta3,tau_epsilon[ix(x,y,z,X,Y,Z)],
							 conc, time, t0[ix(x,y,z,X,Y,Z)],
							 x, y, z, X, Y, Z, T, sigmatheta[ix(x,y,z,X,Y,Z)], aif_settings, spatial, img_mask);

			     // reject kep_new if larger then kep2
			     if(exp(temp) > kep2[ix(x,y,z,X,Y,Z)])
			     {
			    	 temp = log(kep[ix(x,y,z,X,Y,Z)]); // keep old kep
			     }
			     if (temp!=log(kep[ix(x,y,z,X,Y,Z)]))
			     {
			    	 //printf("kep updated: old kep %f, new kep %f\n",kep[ix(x,y,z,X,Y,Z)], exp(temp) );
			    	 acc_theta[ix(x,y,z,X,Y,Z)] = acc_theta[ix(x,y,z,X,Y,Z)] + 1;
			    	 kep[ix(x,y,z,X,Y,Z)] = exp(temp);
			     }


			     // update kep2
				 // (ktrans and ktrans2 / kep and kep2 change role)
				 temp = update_theta(log(kep2[ix(x,y,z,X,Y,Z)]), kep2, ktrans2[ix(x,y,z,X,Y,Z)],kep[ix(x,y,z,X,Y,Z)],ktrans[ix(x,y,z,X,Y,Z)],
											 vp[ix(x,y,z,X,Y,Z)], tau_theta_2comp, tau_theta2_2comp,tau_theta3_2comp,tau_epsilon[ix(x,y,z,X,Y,Z)],
											 conc, time, t0[ix(x,y,z,X,Y,Z)],
											 x, y, z, X, Y, Z, T, sigmatheta_2comp[ix(x,y,z,X,Y,Z)], aif_settings, spatial, img_mask);

				 // reject kep2_new if smaller then kep2
				 if(exp(temp) < kep[ix(x,y,z,X,Y,Z)])
				 {
					 temp = log(kep2[ix(x,y,z,X,Y,Z)]); // keep old kep
				 }

				 // keep kep2_new
				 if (temp!=log(kep2[ix(x,y,z,X,Y,Z)]))
				 {
					 acc_theta_2comp[ix(x,y,z,X,Y,Z)] = acc_theta_2comp[ix(x,y,z,X,Y,Z)] + 1;
					 kep2[ix(x,y,z,X,Y,Z)] = exp(temp);
				 }


			     // update Ktrans (gamma = log(Ktrans))
			     temp = update_gamma(log(ktrans[ix(x,y,z,X,Y,Z)]), ktrans, kep[ix(x,y,z,X,Y,Z)], ktrans2[ix(x,y,z,X,Y,Z)], kep2[ix(x,y,z,X,Y,Z)],
			    		 vp[ix(x,y,z,X,Y,Z)], tau_gamma, tau_gamma2,  tau_gamma3, tau_epsilon[ix(x,y,z,X,Y,Z)], conc, time, t0[ix(x,y,z,X,Y,Z)],
			    		 x, y, z, X, Y, Z, T, sigmagamma[ix(x,y,z,X,Y,Z)], aif_settings, spatial, img_mask);
			     if (temp!=log(ktrans[ix(x,y,z,X,Y,Z)]))
			     {
			    	 acc_gamma[ix(x,y,z,X,Y,Z)]=acc_gamma[ix(x,y,z,X,Y,Z)]+1;
			    	 ktrans[ix(x,y,z,X,Y,Z)]=exp(temp);
			     }

			     // update Ktrans2 (gamma2 = log(Ktrans2))
				 // (ktrans and ktrans2 / kep and kep2 change role)
				 temp = update_gamma3(log(ktrans2[ix(x,y,z,X,Y,Z)]), ktrans2, kep2[ix(x,y,z,X,Y,Z)], ktrans[ix(x,y,z,X,Y,Z)], kep[ix(x,y,z,X,Y,Z)],
						 vp[ix(x,y,z,X,Y,Z)], tau_gamma_2comp, tau_gamma2_2comp, tau_gamma3_2comp, tau_epsilon[ix(x,y,z,X,Y,Z)], conc, time, t0[ix(x,y,z,X,Y,Z)],
						 x, y, z, X, Y, Z, T, sigmagamma_2comp[ix(x,y,z,X,Y,Z)], aif_settings, spatial, img_mask);
				 if (temp!=log(ktrans2[ix(x,y,z,X,Y,Z)]))
				 {
					 acc_gamma_2comp[ix(x,y,z,X,Y,Z)]=acc_gamma_2comp[ix(x,y,z,X,Y,Z)]+1;
					 ktrans2[ix(x,y,z,X,Y,Z)]=exp(temp);
				 }

			     // update vp
			     if (vpupdate==3)
			     {
			    	 temp = update_vp(vp[ix(x,y,z,X,Y,Z)], kep[ix(x,y,z,X,Y,Z)],  ktrans[ix(x,y,z,X,Y,Z)], kep2[ix(x,y,z,X,Y,Z)], ktrans2[ix(x,y,z,X,Y,Z)],
									 tau_vp,  tau_vp2,  tau_epsilon[ix(x,y,z,X,Y,Z)], conc, time,  t0[ix(x,y,z,X,Y,Z)],
									 x, y, z, X, Y, Z, T, sigmavp[ix(x,y,z,X,Y,Z)], aif_settings, a_vp, b_vp);
			    	 if (temp!=vp[ix(x,y,z,X,Y,Z)])
			    	 {
			    		 acc_vp[ix(x,y,z,X,Y,Z)] = acc_vp[ix(x,y,z,X,Y,Z)]+1;
			    		 vp[ix(x,y,z,X,Y,Z)] = temp;
			    	 }
			     }

			     // update of tau_epsilon for voxelwise fit
				 if(spatial==0)
				 {
			    	 double taunew;
			    	 taunew=update_tau_epsilon_voxelwise(a_epsilon, b_epsilon, conc, x, y, z, vp[ix(x,y,z,X,Y,Z)],
									 ktrans[ix(x,y,z,X,Y,Z)], kep[ix(x,y,z,X,Y,Z)], ktrans2[ix(x,y,z,X,Y,Z)], kep2[ix(x,y,z,X,Y,Z)],
									 t0[ix(x,y,z,X,Y,Z)], time, X, Y, Z, T, aif_settings);
			    	 if(taunew >0)
			    	 {
			    		 tau_epsilon[ix(x,y,z,X,Y,Z)]=taunew;
			    	 }

				 }

			 } //endif img_mask[]==1
		 
		 } // end iterate over all voxels
					
		 // update tau_epsilon, if sigma2 assumed same in all voxels
		 // update tau_epsilon (inverse noise variance)
		 if (uptauep>=1)
		 {
			 update_tau_epsilon_all(tau_epsilon, a_epsilon, b_epsilon, conc, img_mask, vp,
										 ktrans,kep,ktrans2,kep2,t0,time, X, Y, Z, T, aif_settings);
		 }

		 // update tau_theta and tau_gamma
		 if (uptauep==2)
		 {
		     update_tau_theta(tau_theta, tau_theta2, tau_theta3, a_theta, a_theta3, b_theta, b_theta3, img_mask, kep, X, Y, Z, N, spatial);

		     update_tau_theta(tau_gamma, tau_gamma2, tau_gamma3, a_gamma, a_gamma3, b_gamma, b_gamma3, img_mask, ktrans, X, Y, Z, N, spatial);


		     update_tau_theta(tau_theta_2comp, tau_theta2_2comp, tau_theta3_2comp, a_theta, a_theta3, b_theta, b_theta3,
			 					 img_mask, kep2, X, Y, Z, N, spatial);

		     update_tau_theta(tau_gamma_2comp, tau_gamma2_2comp, tau_gamma3_2comp, a_gamma, a_gamma3, b_gamma, b_gamma3,
			 				 img_mask, ktrans2, X, Y, Z, N, spatial);
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
			int count7 = 0;
			int count8 = 0;
			int tu=0;
			for (int i=0; i<N; i++)
			{
				if (img_mask[i]==1)
				{
					count2 += acc_theta[i];
					count3 += acc_gamma[i];
					count4 += acc[i];
					count6 += acc_vp[i];
					count7 += acc_theta_2comp[i];
					count8 += acc_gamma_2comp[i];

					if (gemupdate==0)
					{
						tu = space_tune(sigmatheta,(double)acc_theta[i]/(double)acc[i],i);
						tu+= space_tune(sigmagamma,(double)acc_gamma[i]/(double)acc[i],i);

						tu+= space_tune(sigmatheta_2comp,(double)acc_theta_2comp[i]/(double)acc[i],i);
						tu+= space_tune(sigmagamma_2comp,(double)acc_gamma_2comp[i]/(double)acc[i],i);


						if (vpupdate>=2)
						{
							tu+= space_tune(sigmavp,(double)acc_vp[i]/(double)acc[i],i);
						}
					}
					else
					{
						if (nulleins()>0.5)
						{
							tu = space_tune(sigmatheta,(double)acc_theta[i]/(double)acc[i],i);
							tu+= space_tune(sigmatheta_2comp,(double)acc_theta_2comp[i]/(double)acc[i],i);
						}
						else
						{
							tu = space_tune(sigmagamma,(double)acc_gamma[i]/(double)acc[i],i);
							tu+= space_tune(sigmagamma_2comp,(double)acc_gamma_2comp[i]/(double)acc[i],i);
						}
					}
					if (tu>0)
					{
						count5++;
					}

				} //end if (img_mask[i])
			} // end for


			// if the percentage of voxels with tuning is less then space_tunepct
			// -> re-tuning (start next tuning-step, at most respace_tunecycles tuning steps)
			if (100*count5<N1*space_tunepct)
			{
			   respace_tune++;
			   if(!(respace_tune==respace_tunecycles-1))
			   {
				   Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bRe-tuning:\n");
				   iter=0;
			   }
			}
			// else: tune everything again without increasing completed respace_tunecycles
			else
			{
       	   Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%i of ",count5);
           Rprintf("%i done.\n", N1);
			   iter=0;
			}
			for (int i=0; i<N; i++)
			{
			   acc_theta[i]=0;
			   acc_vp[i]=0;
			   acc_gamma[i]=0;
			   acc[i]=0;
			   acc_theta_2comp[i]=0;
			   acc_gamma_2comp[i]=0;

			}

		} // end tuning (if iter==tuning...)

	 } // end for MCMC iterations
	 
	 //Rprintf("End of MCMC iterations\n");

	 PutRNGstate();

     // Calculate deviance
     for (int x=1 ; x<=X  ; x++)
     {
		 for (int y=1  ; y<=Y  ; y++)
		 {
		   for (int z=1; z<=Z; z++)
		   {
			 if(img_mask[ix(x,y,z,X,Y,Z)]==1)
			 {
				 devi[ix(x,y,z,X,Y,Z)] = -2*log_likelihood(tau_epsilon[ix(x,y,z,X,Y,Z)], ktrans[ix(x,y,z,X,Y,Z)], kep[ix(x,y,z,X,Y,Z)],
														   ktrans2[ix(x,y,z,X,Y,Z)], kep2[ix(x,y,z,X,Y,Z)], conc, time,
														   vp[ix(x,y,z,X,Y,Z)], x,y,z,X,Y,Z,T, aif_settings );
			 }
		   }
		 }
     }

     return;
		   
} // end dce_space
