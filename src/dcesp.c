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


void cholesky(int* n, double* matrix,  int* lda, int* info)
{
  F77_CALL(dpotrf)("L", n, matrix, lda, info);
}

void loese(double* A, double* z, int* n, int* lda, int* incx)     
{
  F77_CALL(dtrsv)("L","T","N", n, A, lda, z, incx);
} 

void loese2(double* A, double* z, int* n, int* lda, int* incx)     
{
 F77_CALL(dtrsv)("L","N","N", n, A, lda, z, incx);
} 

void copy(double* A,double* B,int n1,int n2)
{
  int i;
  for (i=0;i<n1*n2;i++)
    {
      A[i]=B[i];
      
    }
}
void transpose(double* A,double* B,int n1,int n2)
{
  int i,j;
  for (i=0;i<n1;i++)
    {
      for (j=0;j<n2;j++)
	{
	  A[j+i*n2]=B[i+j*n1];
	}
    }
}

void multiply(double* E,double* A,double* B,int* n1,int* n2,int* n3)
{
  int i,j,k;

  for (i=0;i<*n1;i++)
    {
      for (j=0;j<*n3;j++)
	{
	 
	    E[i+*n1*j]=0.0;
	    
	  for (k=0;k<*n2;k++)
	    {
	      E[i+*n1*j]=E[i+*n1*j]+A[i+*n1*k]*B[k+*n2*j];
	    }
	  
	}
    }
}

void multiplyvA(double* E,double* A,double* B, int* n2, int* n3)
{
  int j,k;
  

  for (j=0;j<*n3;j++)
    {
      //Rprintf("%i ",j);
      
	E[j]=0.0;
     
      for (k=0;k<*n2;k++)
	{
	  E[j]=E[j]+A[k]*B[k+*n2*j];
	}
     
    }
}

void add(double* E,double* A,double* B,int n1,int n2)
{
  int i;
  for (i=0;i<n1*n2;i++)
    {
	  E[i]=A[i]+B[i];
    }
}

void R(double* R, double* v, int rw, int n){
  int i,j;

  for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
	{
	  R[i+n*j]=0.0;
	}
    }
  if (rw==0)
    {
      for (i=0;i<(n-rw);i++)
	{
	  R[i+n*i]=v[i];
	}
    }
  if (rw==1)
    {
      for (i=0;i<(n-rw);i++)
	{
	  R[i+n*i]=R[i+n*i]+v[i];
	  R[i+1+n*(i+1)]=R[i+1+n*(i+1)]+v[i];
	  R[i+1+n*i]=R[i+1+n*i]-v[i];
	  R[i+n*(i+1)]=R[i+n*(i+1)]-v[i];
	}
    }
  if (rw==2)
    {
      for (i=0;i<(n-rw);i++)
	{
	  R[i+n*i]=R[i+n*i]+v[i];
	  R[i+1+n*(i+1)]=R[i+1+n*(i+1)]+4*v[i];
	  R[i+2+n*(i+2)]=R[i+2+n*(i+2)]+v[i];
	  R[i+1+n*i]=R[i+1+n*i]-2*v[i];
	  R[i+n*(i+1)]=R[i+n*(i+1)]-2*v[i];
	  R[i+2+n*i]=R[i+2+n*i]+v[i];
	  R[i+n*(i+2)]=R[i+n*(i+2)]+v[i];
	  R[i+2+n*(i+1)]=R[i+2+n*(i+1)]-2*v[i];
	  R[i+1+n*(i+2)]=R[i+1+n*(i+2)]-2*v[i];
	}
    }

}

void Rv(double* R, double v, int rw, int n){
  int i,j;

  for (i=0;i<n;i++)
    {
      for (j=0;j<n;j++)
	{
	  R[i+n*j]=0.0;
	}
    }
  if (rw==0)
    {
      for (i=0;i<(n-rw);i++)
	{
	  R[i+n*i]=v;
	}
    }
  if (rw==1)
    {
      for (i=0;i<(n-rw);i++)
	{
	  R[i+n*i]=R[i+n*i]+v;
	  R[i+1+n*(i+1)]=R[i+1+n*(i+1)]+v;
	  R[i+1+n*i]=R[i+1+n*i]-v;
	  R[i+n*(i+1)]=R[i+n*(i+1)]-v;
	}
    }
  if (rw==2)
    {
      for (i=0;i<(n-rw);i++)
	{
	  R[i+n*i]=R[i+n*i]+v;
	  R[i+1+n*(i+1)]=R[i+1+n*(i+1)]+4*v;
	  R[i+2+n*(i+2)]=R[i+2+n*(i+2)]+v;
	  R[i+1+n*i]=R[i+1+n*i]-2*v;
	  R[i+n*(i+1)]=R[i+n*(i+1)]-2*v;
	  R[i+2+n*i]=R[i+2+n*i]+v;
	  R[i+n*(i+2)]=R[i+n*(i+2)]+v;
	  R[i+2+n*(i+1)]=R[i+2+n*(i+1)]-2*v;
	  R[i+1+n*(i+2)]=R[i+1+n*(i+2)]-2*v;
	}
    }

}
int ix5(int a,int b,int c,int d,int e,int A,int B,int C,int D,int E)
{
  return(e*A*B*C*D+d*A*B*C+c*A*B+b*A+a);
}
int ix4(int a,int b,int c,int d, int A,int B,int C,int D)
{
  return(d*A*B*C+c*A*B+b*A+a);
}
int ix3(int a,int b,int c,int A,int B,int C)
{
  return(c*A*B+b*A+a);
}
int ix2(int a,int b, int A,int B)
{
  return(b*A+a);
}
	
void vector5(double* res,double* arra,int l,int n,int a,int b,int c,int d,
	     int A, int B, int C, int D, int E){
  int i;
  if (n==1)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix5(i,a,b,c,d,A,B,C,D,E)];
	}
    }
  if (n==2)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix5(a,i,b,c,d,A,B,C,D,E)];
	}
    }
  if (n==3)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix5(a,b,i,c,d,A,B,C,D,E)];
	}
    }
  if (n==4)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix5(a,b,c,i,d,A,B,C,D,E)];
	}
    }
  if (n==5)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix5(a,b,c,d,i,A,B,C,D,E)];
	}
    }

}

void vector4(double* res,double* arra,int l,int n,int a,int b,int c, int A, int B, int C, int D){
  int i;
  if (n==1)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix4(i,a,b,c,A,B,C,D)];
	}
    }
  if (n==2)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix4(a,i,b,c,A,B,C,D)];
	}
    }
  if (n==3)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix4(a,b,i,c,A,B,C,D)];
	}
    }
  if (n==4)
    {
      for (i=0;i<l;i++)
	{
	  res[i]=arra[ix4(a,b,c,i,A,B,C,D)];
	}
    }
}

void draw_normal(double* x, double* L, int* n, double* temp)
{

  int i,c,d;
  c=1;
  d=0;
  cholesky(n, L, n, &d);
  loese2(L, x, n, n ,&c);
  loese(L, x, n, n ,&c);
  gausssample(temp,n);
  loese(L,temp,n, n ,&c);
  for (i=0; i<*n; i++)
  {
    x[i]=x[i]+temp[i];
  }

}

void dce_spline_run(int* NRI, int* thin, int* dims, double* data,
	 double* tau, double* tauepsilon, double* D,
	 int* rw, double* beta, double* hyper, int* p,
  double* tmpvec,  double* tmpvec1,
  double* tmpvec2,  double* Mx,
  double* Mx1,  double* Mx2,
  double* Mx3,	 double* Dt)
{
  GetRNGstate();
  int i,j,x,y,z,k;
  double tmp,tmp1;
  int X=dims[0];
  int Y=dims[1];
  int Z=dims[2];
  int T=dims[3];


  transpose(Dt,D,T,*p);
 
  for (i=0; i<*NRI; i++)
    {
/*       if (fmod(i,50)==0){ */
/* 	Rprintf("Iteration %i\n",i+1); */
/*       } */
      for (j=0; j<*thin; j++)
	{
	  for (x=0; x<dims[0]; x++)
	    {
	      for (y=0; y<dims[1]; y++)
		{
		  for (z=0; z<dims[2]; z++)
		    {
		      
		      vector5(tmpvec1,tau,*p-*rw,4,x,y,z,i,X,Y,Z,*p-*rw,*NRI);
		      //vector5(tmpvec,tauepsilon,T,4,x,y,z,i,X,Y,Z,T,*NRI);
		      tmp=tauepsilon[ix4(x,y,z,i,X,Y,Z,*NRI)];
		      Rv(Mx,tmp,0,T);
		      multiply(Mx2,Dt,Mx,p,&T,&T);
		      multiply(Mx3,Mx2,D,p,&T,p);
		      R(Mx1,tmpvec1,*rw,*p);
		      add(Mx2,Mx3,Mx1,*p,*p); //Mx2 is precision matrix
		      //vector5(tmpvec,tauepsilon,T,4,x,y,z,i,X,Y,Z,T,*NRI);
		      //R(Mx,tmpvec,0,T);
		      Rv(Mx,tmp,0,T);
		      vector4(tmpvec,data,T,4,x,y,z,X,Y,Z,T);
		      multiplyvA(tmpvec1,tmpvec,Mx,&T,&T);
		      multiplyvA(tmpvec,tmpvec1,D,&T,p); //tmpvec is canonical mean
		      draw_normal(tmpvec,Mx2,p,tmpvec2);
		    
		      for (k=0; k<*p; k++)
			{
			  beta[ix5(x,y,z,k,i,X,Y,Z,*p,*NRI)]=tmpvec[k];
			}
		      
		      if (*rw==1)
			{
			  tmpvec1[0]=hyper[0]+0.5;
			  tmpvec1[*p-2]=hyper[0]+0.5;
			  for (k=1;k<(*p-*rw-1);k++)
			    {
			      tmpvec1[k]=hyper[0]+1;
			    }
			  for (k=0;k<(*p-1);k++)
			    {
			      tmpvec2[k]=hyper[1]+0.5*(beta[ix5(x,y,z,k,i,X,Y,Z,*p,*NRI)]-beta[ix5(x,y,z,k+1,i,X,Y,Z,*p,*NRI)])*(beta[ix5(x,y,z,k,i,X,Y,Z,*p,*NRI)]-beta[ix5(x,y,z,k+1,i,X,Y,Z,*p,*NRI)]);
			    }
			}
		      if (*rw==2)
			{
			  tmpvec1[0]=hyper[0]+0.5;
			  tmpvec1[1]=hyper[0]+1;
			  tmpvec1[*p-3]=hyper[0]+0.5;
			  tmpvec1[*p-4]=hyper[0]+1;
			  for (k=2;k<(*p-4);k++)
			    {
			      tmpvec1[k]=hyper[0]+3/2;
			    }
			  for (k=0;k<(*p-2);k++)
			    {
			      tmp=beta[ix5(x,y,z,k,i,X,Y,Z,*p,*NRI)]-2*beta[ix5(x,y,z,k+1,i,X,Y,Z,*p,*NRI)]+beta[ix5(x,y,z,k+2,i,X,Y,Z,*p,*NRI)];
			      tmpvec2[k]=hyper[1]+0.5*tmp*tmp;
			    }
			}
		      for (k=0;k<(*p-*rw);k++)
			{
			  tau[ix5(x,y,z,k,i,X,Y,Z,*p-*rw,*NRI)]=RNDGAM(tmpvec1[k],tmpvec2[k]);
			}
		      
		      vector5(tmpvec,beta,*p,4,x,y,z,i,X,Y,Z,*p,*NRI); 
		      multiplyvA(tmpvec1,tmpvec,Dt,p,&T);
		      tmp1=0.0; 
		      for (k=0;k<T;k++)
			{ 
			  tmp=data[ix4(x,y,z,k,X,Y,Z,T)]-tmpvec1[k];
			  tmp1+=tmp*tmp;
			} 
  		      tauepsilon[ix4(x,y,z,i,X,Y,Z,*NRI)]=RNDGAM(hyper[2]+(double)T/2.0,hyper[3]+0.5*tmp1);
	      
		    }
		}
	    }
	}
      
      if ((i+1)!=*NRI)
	{
	  for (x=0; x<dims[0]; x++)
	    {
	      for (y=0; y<dims[1]; y++)
		{
		  for (z=0; z<dims[2]; z++)
		    {
		      for (k=0; k<(*p-*rw);k++)
			{
			  tau[ix5(x,y,z,k,i+1,X,Y,Z,*p-*rw,*NRI)]=tau[ix5(x,y,z,k,i,X,Y,Z,*p-*rw,*NRI)];
			}
		      for (k=0; k<*p;k++)
			{
			  beta[ix5(x,y,z,k,i+1,X,Y,Z,*p,*NRI)]=beta[ix5(x,y,z,k,i,X,Y,Z,*p,*NRI)];
			}
		      tauepsilon[ix4(x,y,z,i+1,X,Y,Z,*NRI)]=tauepsilon[ix4(x,y,z,i,X,Y,Z,*NRI)];

		    }
		}
	    }
	}
    
    }
  PutRNGstate();

}

	
