#include <iostream>
#include <fstream>
#include <complex>
#include <iomanip>
#include <time.h>
#include <stdlib.h>


#define complex complex<double>
using namespace std;

#define do_once  {static int flag_once=0; if (flag_once==0) {flag_once=1;
#define end_do_once ;};}


#define for1(i, N) for (int i=0; i<N; i++) 
#define for2(i, j, N) for (int i=0; i<N; i++) for (int j=0; j<N; j++)
#define for3(i, j, k, N) for (int i=0; i<N; i++) for (int j=0; j<N; j++) for (int k=0; k<N; k++)

int sgn(double x) {if (x>0) return 1; if (x<0) return -1; return 0;}
double sqr(double x) {return x*x;}
complex sqr(const complex & x) {return x*x;}

double norm2(const complex & x) {return real(x)*real(x)+imag(x)*imag(x);}


void swing(int & x, int & y) {int a=x; x=y; y=a;}
void swing(double & x, double & y) {double  a=x; x=y; y=a;}
void swing(complex & x, complex & y) {complex a=x; x=y; y=a;}
void swing(double * & x, double * & y) {double * a=x; x=y; y=a;}
void swing(complex * & x, complex * & y) {complex * a=x; x=y; y=a;}
 

int min(int & x, int &y) {if (x<y) return x; else return y;}
double min(double x, double y) {if (x<y) return x; else return y;}

int max(int & x, int &y) {if (x>y) return x; else return y;}

complex I(0,1); double Pi=4*atan(1);

int initial_random=time(NULL);
int INT_RANDOM=initial_random;//  //from 1 to 2^31-1

int int_rnd()
{
	 int k=INT_RANDOM/127773;
	 INT_RANDOM=16807*(INT_RANDOM-k*127773)-2836*k;
	 if (INT_RANDOM<0) INT_RANDOM+=2147483647;
	 return INT_RANDOM;
;}

int INT_RANDOM_AUX_VARIABLE=int_rnd()+int_rnd()+int_rnd()+int_rnd();  //randomising...

double rnd ()
{
	return int_rnd()/2147483647.0
;}

int rnd (int k)
{
	int d=2147483647%k, d1=2147483647-d, r=int_rnd();
	if (r>=d1) return rnd(k);
	return r/(d1/k)
;}

complex rnd_gauss2()
{
  double R=sqrt(-log(rnd())), phi=2*Pi*rnd(); return R*exp(I*phi);  
//Gauss::   double R=sqrt(-2.*log(rnd())), phi=2*Pi*rnd(); return R*exp(I*phi);  
}

int ** new_int2(int n1, int n2)
{
  int ** r=new int *[n1];
  for (int i=0; i<n1; i++)  r[i]=new int [n2];   
  return r;
}

void delete_int2(int ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}


double ** new_double2(int n1, int n2)
{
  double ** r=new double *[n1];
  for (int i=0; i<n1; i++)  r[i]=new double [n2];   
  return r;
}

void delete_double2(double ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}

double ** new_ldouble2(int n1, int n2)
{
  double ** r=new double *[n1];
  for (int i=0; i<n1; i++)  r[i]=new double [n2];   
  return r;
}

void delete_ldouble2(double ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}


complex ** new_complex2(int n1, int n2)
{
  complex ** r=new complex *[n1];
  for (int i=0; i<n1; i++)  r[i]=new complex [n2];   
  return r;
}

void delete_complex2(complex ** &r, int n1, int n2)
{
  for (int i=0; i<n1; i++)  delete [] r[i];   
  delete [] r;
}


int ***  new_int3(int n1, int n2, int n3)
{
  int *** r;
  r=new int **[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new int *[n2];
    for (int j=0; j<n2; j++) r[i][j]=new int [n3];
  }
  return r;
}

void delete_int3(int *** &r, int n1, int n2, int n3)
{
  for (int i=0; i<n1; i++)  
  {
    for (int j=0; j<n2; j++) delete [] r[i][j];
    delete [] r[i];  
  }
  delete [] r;
}

double *** new_double3(int n1, int n2, int n3)
{
  double *** r=new double **[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new double *[n2];
    for (int j=0; j<n2; j++) r[i][j]=new double [n3];
  }
  return r;
}

void delete_double3(double *** &r, int n1, int n2, int n3)
{
  for (int i=0; i<n1; i++)  
  {
    for (int j=0; j<n2; j++) delete [] r[i][j];
    delete [] r[i];  
  }
  delete [] r;
}

complex *** new_complex3(int n1, int n2, int n3)
{
  complex *** r=new complex **[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new complex *[n2];
    for (int j=0; j<n2; j++) r[i][j]=new complex [n3];
  }
  return r;
}

void delete_complex3(complex *** &r, int n1, int n2, int n3)
{
  for (int i=0; i<n1; i++)  
  {
    for (int j=0; j<n2; j++) delete [] r[i][j];
    delete [] r[i];  
  }
  delete [] r;
}


int **** new_int4(int n1, int n2, int n3, int n4)
{
  int **** r=new int ***[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new int **[n2];
    for (int j=0; j<n2; j++) 
    {
      r[i][j]=new int * [n3];
      for (int l=0; l<n3; l++) r[i][j][l]=new int [n4];
    }
  }
  return r;
}


complex **** new_complex4(int n1, int n2, int n3, int n4)
{
  complex **** r=new complex ***[n1];
  for (int i=0; i<n1; i++) 
  {
    r[i]=new complex **[n2];
    for (int j=0; j<n2; j++) 
    {
      r[i][j]=new complex * [n3];
      for (int l=0; l<n3; l++) r[i][j][l]=new complex [n4];
    }
  }
  return r;
}



void delete_complex4(complex **** &r, int n1, int n2, int n3, int n4)
{
  for (int i=0; i<n1; i++)  
  {
    for (int j=0; j<n2; j++) 
    {
        for(int l=0; l<n3; l++) delete [] r[i][j][l];
        delete [] r[i][j];
    }
    delete [] r[i];  
  }
  delete [] r;
}

void Inverse(complex ** &a, int size) //Gauss with partial pivoting
{
	complex ** r=new complex * [size];
	{
		for (int i=0;i<size;i++)
		{
			r[i]=new complex [size]; 
			for (int j=0; j<size; j++) r[i][j]=0.; r[i][i]=1.;
		;}
	;}
	
	for (int i=0; i<size; i++)
	{
		//pivoting
      {
		double fmax=-1.; int jmax;
		for (int j=i; j<size; j++)
			if ( fmax<norm2(a[j][i]) ) {jmax=j; fmax=norm2(a[j][i]);}
		complex * aux; 
		aux=a[i]; a[i]=a[jmax]; a[jmax]=aux;
		aux=r[i]; r[i]=r[jmax]; r[jmax]=aux;
		if (norm2(a[i][i])<=0.) {cout<<"LInverse!"; return;}
      ;}
		
		//main body
	
#pragma omp parallel for
		for (int j=0; j<i; j++)
		{
			complex f=a[j][i]/a[i][i];a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f; 
			for (int l=0; l<size; l++) r[j][l]-=r[i][l]*f;
			}
		;}
		{
			complex f=1./a[i][i];
			{for (int l=i; l<size; l++) a[i][l]*=f;}
			{for (int l=0; l<size; l++) r[i][l]*=f;}
		;}
#pragma omp parallel for
		for (int j=i+1; j<size; j++)
		{
			complex f=a[j][i]/a[i][i]; a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f;
			for (int l=0; l<size; l++) r[j][l]-=r[i][l]*f;
			}
		;}
		
        ;}
	
	{for (int i=0; i<size; i++) for (int j=0; j<size; j++) a[i][j]=r[i][j];}
	{for (int i=0; i<size; i++) delete [] r[i]; delete [] r;}
;}


void div_left(complex ** a, complex *  r, int size) //a=r^{-1}*a, Gauss with partial pivoting
//destroys data in a
{

    for (int i=0; i<size; i++)
	{
		//pivoting
      {
		double fmax=-1.; int jmax=i;
		for (int j=i; j<size; j++)
			if ( fmax<norm2(a[j][i]) ) {jmax=j; fmax=norm2(a[j][i]);}
		swing(a[i], a[jmax]); swing(r[i], r[jmax]);
		if (norm2(a[i][i])<=0.) {cout<<"Error in div_left "; return;}
      ;}
		
		//main body
	
#pragma omp parallel for
		for (int j=0; j<i; j++)
		{
			complex f=a[j][i]/a[i][i];a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f; 
			r[j]-=r[i]*f;
			}
		;}
		{
			complex f=1./a[i][i];
			{for (int l=i; l<size; l++) a[i][l]*=f;}
			r[i] *=f;
		;}
#pragma omp parallel for
		for (int j=i+1; j<size; j++)
		{
			complex f=a[j][i]/a[i][i]; a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f;
			r[j]-=r[i]*f;
			}
		;}
    }	
    
}



bool div_left(double ** a, double *  r, int size) //a=r^{-1}*a, Gauss with partial pivoting
//destroys data in a
{

    for (int i=0; i<size; i++)
	{
		//pivoting
      {
		double fmax=-1.; int jmax=i;
		for (int j=i; j<size; j++)
			if ( fmax<abs(a[j][i]) ) {jmax=j; fmax=abs(a[j][i]);}
		swing(a[i], a[jmax]);     swing(r[i], r[jmax]);  
		if (abs(a[i][i])<=0.) return false;//{cout<<"Error in div_left "; return;}
      ;}
		//main body
	
#pragma omp parallel for
		for (int j=0; j<i; j++)
		{
			double f=a[j][i]/a[i][i];a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f; 
			r[j]-=r[i]*f;
			}
		;}
		{
			double f=1./a[i][i];
			{for (int l=i; l<size; l++) a[i][l]*=f;}
			r[i] *=f;
		;}
#pragma omp parallel for
		for (int j=i+1; j<size; j++)
		{
			double f=a[j][i]/a[i][i]; a[j][i]=0;
			{
			for (int l=i+1; l<size; l++) a[j][l]-=a[i][l]*f;
			r[j]-=r[i]*f;
			}
		;}
    }	
    return true;    
}



bool div_conj_grad(double ** A, double *  b, double *x, int size)//could be optimized
{
    static int size_old=0; static double *p=NULL, * r=NULL, *Ap=NULL;
    if (size>size_old) {delete [] p; delete [] r; delete [] Ap; p=new double  [size]; r=new double  [size]; Ap=new double [size]; size_old=size;}
    
     
    for1(j, size) {x[j]=0;r[j]=b[j];p[j]=b[j];}
    
    for1(i,5*size)
    {  
        double r2=0; for1(j, size) r2+=r[j]*r[j];
        for1(j, size) {Ap[j]=0; for1(j_, size) Ap[j]+=A[j][j_]*p[j_];}
        
        double p2=0; for1(j,size) p2+=p[j]*Ap[j];
        double a=r2/p2;  
        for1(j, size) x[j]+=a*p[j];  

        for1(j, size) r[j]-= a*Ap[j];
        double r2_=0; for1(j, size) r2_+=r[j]*r[j];
        cout<<i<<"  "<<r2<<"  ";
        if (r2_<1e-20) return true;
        double b=r2_/r2;
        for1(j, size) p[j]=r[j]+b*p[j];
    }
    return false;
}


void EigenJacobi (double ** h, double ** a, double * e, int n, double accuracy=1e-10, int max_sweep=100)
{
  int flag; double ** v=new double * [n]; {for (int i=0; i<n; i++) v[i]=new double [n];}
  
  {for (int i=0; i<n; i++) for (int j=0; j<n; j++) {a[i][j]=h[i][j]; if (i==j) v[i][j]=1; else v[i][j]=0;};}

	do
	{
		max_sweep--;
		flag=0;
		for (int pp=n-1; pp>=1; pp--)
		for (int q=0; q<n-pp; q++)
		{
			int p=pp+q;
			if (abs(a[p][q])>accuracy)
			{
				flag=1;
				double theta=(a[q][q]-a[p][p])/(2*a[p][q]);
				double t=1/(abs(theta)+sqrt(theta*theta+1));if (theta<0) t=-t;
				double c=1/sqrt(t*t+1), s=t*c, tau=s/(1+c);
				for (int r=0; r<n; r++)
				{
					if ((r!=p) && (r!=q))
					{
						double a1=a[r][p]-s*(a[r][q]+tau*a[r][p]);
						double a2=a[r][q]+s*(a[r][p]-tau*a[r][q]);
						a[r][p]=a1; a[p][r]=a1;
						a[r][q]=a2; a[q][r]=a2;


					;}
						double v1=c*v[p][r]-s*v[q][r];
						double v2=s*v[p][r]+c*v[q][r];
						v[p][r]=v1; v[q][r]=v2;
				;}
				double ap=a[p][p]-t*a[p][q], aq=a[q][q]+t*a[p][q];
				a[p][p]=ap; a[q][q]=aq;
				a[p][q]=0;  a[q][p]=0;
			;}
		;}
	;}
	while( (flag>0) && (max_sweep>0) );
	
	
	for (int i=0; i<n; i++) e[i]=a[i][i];

/*
	for (int i=0; i<n-1; i++)
	for (int j=i+1; j<n; j++)
		if (e[j]<e[i])
		{
			double a=e[j]; e[j]=e[i]; e[i]=a;
			double * b=v[j]; v[j]=v[i];v[i]=b;
		;}
*/
	{for (int i=0; i<n; i++) for (int j=0; j<n; j++) a[i][j]=v[i][j];}
	


	for (int i=0; i<n; i++) //normalization
	{
		double s=0;
		for (int j=0; j<n; j++) s+=a[i][j]*a[i][j];
		s=sqrt(s);
		if (abs(s)>accuracy*accuracy) {for (int j=0; j<n; j++) a[i][j]/=s;}
		else {for (int j=0; j<n; j++) a[i][j]=1/sqrt(n);}
	;}
	
	{for (int i=0; i<n; i++) delete [] v[i]; delete [] v;}
;}


void EigenBlock (double ** h, double ** a, double * e, int n, int * block, int n_block, double accuracy=1e-10, int max_sweep=100)
{
    for2(j1,j2,n) a[j1][j2]=0;
    int * ib=new int [n];
    double ** hb=new_double2(n,n), ** ab=new_double2(n,n), * eb=new double [n];
    for1(b, n_block)
    {
        int nb=0; for1(j, n) if (block[j]==b) {ib[nb]=j; nb++;}
        if (nb!=0)
        {
            for2(j1,j2,nb) hb[j1][j2]=h[ib[j1]][ib[j2]];
            EigenJacobi(hb, ab,eb, nb, accuracy, max_sweep);
            for2(j1,j2,nb) a[ib[j1]][ib[j2]]=ab[j1][j2]; for1(j, nb) e[ib[j]]=eb[j]; 
        }
    }
    delete [] ib; delete_double2(hb,n,n); delete_double2(ab,n,n); delete [] eb;
}



int Eigen_leading (double ** h, double ** a, double * e, int n, double emin=1e-10)
{
    static double * psi=NULL, * psi_=NULL; static int nmax=0;
    if (n>nmax)
    {delete [] psi; delete [] psi_; nmax=n; psi=new double [n]; psi_=new double [n];}
    
    for1(n0, n)
    {
        for1(j, n) psi[j]=rnd()-.5; double s;
        for1(l, 10000)
        {
            for1(j, n) {psi_[j]=0; for1(j1, n) psi_[j]+=h[j][j1]*psi[j1];}
            for1(j, n0){double x=0; for1(j1, n) x+=a[j][j1]*psi_[j1]; for1 (j1, n) psi_[j1]-=x*a[j][j1];}
            s=0; for1(j, n) s+=psi_[j]*psi_[j]; s=sqrt(s); for1(j, n) psi_[j]/=s;
            double r=0, r_=0; for1(j, n) {r+=sqr(psi[j]-psi_[j]); r_+=sqr(psi[j]+psi_[j]); psi[j]=psi_[j];} 
            if (l>990) {for1(j, n) cout<<psi[j]<<"  "; cout<<"\n";}
            
            if (r_<r) s=-s;
            
            if (r<1e-20 || r_<1e-20) {   break;}                                                       
                                                       
        }
        for1(j, n) a[n0][j]=psi[j]; e[n0]=s;
        if (abs(s)<emin) {return  n0+1;}
        
    }
    
    return n;
    
    
}
                  
                  
int Eigen_leading (complex ** h, complex ** a, double * e, int n, double emin=1e-10)
{
    static complex * psi=NULL, * psi_=NULL; static int nmax=0;
    if (n>nmax)
    {delete [] psi; delete [] psi_; nmax=n; psi=new complex [n]; psi_=new complex [n];}
    
    for1(n0, n)
    {
        for1(j, n) psi[j]=rnd()-.5+I*(rnd()-.5); double s;
        for1(l, 10000)
        {
            for1(j, n) {psi_[j]=0; for1(j1, n) psi_[j]+=h[j][j1]*psi[j1];}
            for1(j, n0){complex x=0; for1(j1, n) x+=conj(a[j][j1])*psi_[j1]; for1 (j1, n) psi_[j1]-=x*a[j][j1];}
            s=0; for1(j, n) s+=norm2(psi_[j]); s=sqrt(s); for1(j, n) psi_[j]/=s;
            
            double r=0, r_=0; for1(j, n) {r+=norm2(psi[j]-psi_[j]); r_+=norm2(psi[j]+psi_[j]); psi[j]=psi_[j];} 
            if (l>9990) {for1(j, n) cout<<psi[j]<<"  "; cout<<"\n";}
            
            if (r_<r) s=-s;
            
            if (r<1e-18 || r_<1e-18) {   break;}                                                       
                                                       
        }
        for1(j, n) a[n0][j]=psi[j]; e[n0]=s;
        if (abs(s)<emin) return n0+1;
        
    }
    
    return n;
    
    
}
                  
                  
