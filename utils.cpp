//
//  utils.cpp
//  Project1_555
//
//  Created by Bing Hong Fu on 2/2/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

double ran1(int *idum) {
    static long ix1,ix2,ix3;
    static double r[98];
    double temp;
    static int iff=0;
    int j;
    void nrerror(char *);
    
    if (*idum < 0 || iff == 0) {
        iff=1;
        ix1=(IC1-(*idum)) % M1;
        ix1=(IA1*ix1+IC1) % M1;
        ix2=ix1 % M2;
        ix1=(IA1*ix1+IC1) % M1;
        ix3=ix1 % M3;
        for (j=1;j<=97;j++) {
            ix1=(IA1*ix1+IC1) % M1;
            ix2=(IA2*ix2+IC2) % M2;
            r[j]=(ix1+ix2*RM2)*RM1;
        }
        *idum=1;
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;
    j=1 + ((97*ix3)/M3);
    if (j > 97 || j < 1) nrerror("RAN1: This cannot happen.");
    temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;
    return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3




#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
/*Long period (> 2 x 10^18 ) random number generator of L'Ecuyer with Bays-Durham shuffle
 and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
 the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
 idum between successive deviates in a sequence. RNMX should approximate the largest floating
 value that is less than 1.*/
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if (*idum <= 0) {               //Initialize.
        if (-(*idum) < 1) *idum=1;  //Be sure to prevent idum = 0.
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {   //Load the shu.e table (after 8 warm-ups).
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    
    k=(*idum)/IQ1;                  //Start here when not initializing.
    *idum=IA1*(*idum-k*IQ1)-k*IR1;  //Compute idum=(IA1*idum) % IM1 without overflows by Schrage's method.
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;  //Compute idum2=(IA2*idum) % IM2 likewise.
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;                      //Will be in the range 0..NTAB-1.
    iy=iv[j]-idum2;                 //Here idum is shu.ed, idum and idum2 are combined to generate output.
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX; //Because users don't expect endpoint values.
    else return temp;
}


#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX




int uni_int(int *idum,int m) {
    double u;
    while ((u=ran1(idum))>=1); return (int) floor(m*u);
}

int uni_int2(long int *idum,int m) {
    double u;
    while ((u=ran2(idum))>=1); return (int) floor(m*u);
}


double gauss1(int *idum) {
    static int iset=0;
    static double gset;
    double fac,r,v1,v2;
    //double ran1(int *);
    
    if  (iset == 0) {
        do {
            v1=2.0*ran1(idum)-1.0;
            v2=2.0*ran1(idum)-1.0;
            r=v1*v1+v2*v2;
        } while (r >= 1.0);
        fac=sqrt(-2.0*log(r)/r);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}



double gauss2(long int *idum) {
    static int iset=0;
    static double gset;
    double fac,r,v1,v2;
    
    if  (iset == 0) {
        do {
            v1=2.0*ran2(idum)-1.0;
            v2=2.0*ran2(idum)-1.0;
            r=v1*v1+v2*v2;
        } while (r >= 1.0);
        fac=sqrt(-2.0*log(r)/r);
        gset=v1*fac;
        iset=1;
        return v2*fac;
    } else {
        iset=0;
        return gset;
    }
}


void nrerror(char *errortext)
{
    fprintf(stderr,"Run-time error has occurred ...\n");
    fprintf(stderr,"%s\n",errortext);
    fprintf(stderr,"... now exiting to the system... \n");
    exit(1);
}


double *vector(int ilow,int ihigh)
{
    double *v;
    v=(double *)malloc((size_t) ((ihigh-ilow+1)*sizeof(double)));
    if (!v) nrerror("allocation failure in vector()");
    return v-ilow;
}


void free_vector(double *v,int ilow,int ihigh)
{
    free((void*) (v+ilow));
}

double **matrix(int ilow,int ihigh,int jlow,int jhigh)
{
    int i;
    double **m;
    
    /* Allocate pointers to the rows  */
    m=(double **)malloc((size_t) (ihigh-ilow+1)*sizeof(double*));
    if (!m) nrerror("allocation error in matrix()");
    m -= ilow;
    
    /* allocate space for the rows and set the pointers to them  */
    for (i=ilow;i<=ihigh;i++)
    {
        m[i]=(double *)malloc((size_t) (jhigh-jlow+1)*sizeof(double));
        if(!m[i]) nrerror("allocation error 2 in dmatrix()");
        m[i] -= jlow;
    }
    return m;
}

void free_matrix(double **m,int ilow,int ihigh,int jlow,int jhigh)
{
    int i;
    
    for(i=ihigh;i>=ilow;i--) free((void*) (m[i]+jlow));
    free((void*) (m+ilow));
}

int *ivector(int ilow,int ihigh)
{
    int *v;
    
    v=(int *)malloc((size_t) ((ihigh-ilow+1)*sizeof(int)));
    if (!v) nrerror("allocation failure in vector()");
    return v-ilow;
}
void free_ivector(int *v,int ilow,int ihigh)
{
    free((void*) (v+ilow));
}

int **imatrix(int ilow,int ihigh,int jlow,int jhigh)
{
    int i;
    int **m;
    
    /* Allocate pointers to the rows  */
    m=(int **)malloc((size_t) (ihigh-ilow+1)*sizeof(int*));
    if (!m) nrerror("allocation error in matrix()");
    m -= ilow;
    
    /* allocate space for the rows and set the pointers to them  */
    for (i=ilow;i<=ihigh;i++)
    {
        m[i]=(int *)malloc((size_t) (jhigh-jlow+1)*sizeof(int));
        if(!m[i]) nrerror("allocation error 2 in dmatrix()");
        m[i] -= jlow;
    }
    return m;
}

void free_imatrix(int **m,int ilow,int ihigh,int jlow,int jhigh)
{
    int i;
    
    for(i=ihigh;i>=ilow;i--) free((void*) (m[i]+jlow));
    free((void*) (m+ilow));
}




char **charmatrix(int ilow,int ihigh,int jlow,int jhigh)
{
    int i;
    char **m;
    
    /* Allocate pointers to the rows  */
    m=(char **)malloc((size_t) (ihigh-ilow+1)*sizeof(char*));
    if (!m) nrerror("allocation error in charmatrix()");
    m -= ilow;
    
    /* allocate space for the rows and set the pointers to them  */
    for (i=ilow;i<=ihigh;i++)
    {
        m[i]=(char *)malloc((size_t) (jhigh-jlow+1)*sizeof(char));
        if(!m[i]) nrerror("allocation error 2 in charmatrix()");
        m[i] -= jlow;
    }
    return m;
}

void free_charmatrix(char **m,int ilow,int ihigh,int jlow,int jhigh)
{
    int i;
    
    for(i=ihigh;i>=ilow;i--) free((void*) (m[i]+jlow));
    free((void*) (m+ilow));
}



double mag_sq(double *v,int ilow,int ihigh)
{
    int i;
    double sum=0;
    for(i=ilow;i<=ihigh;++i) sum+=pow(v[i],2.0);
    return sum;
}

int ipow(int m, int p)
/* do the pow() function for integers                           */
{
    int i,prod;
    if(!p) return 1;
    if(p<0) nrerror("ipow called with a negative exponent");
    prod=m;
    for(i=2;i<=p;++i) prod*=m;
    return prod;
}


double ***tensor(int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh)
{
    int i,j;
    double ***m;
    
    /* Allocate pointers to the matrices  */
    m=(double ***)malloc((size_t) (ihigh-ilow+1)*sizeof(double**));
    if (!m) nrerror("allocation error 1 in tensor()");
    m -= ilow;
    
    /* allocate pointers to the rows */
    for (i=ilow;i<=ihigh;i++)
    {
        m[i]=(double **)malloc((size_t) (jhigh-jlow+1)*sizeof(double*));
        if(!m[i]) nrerror("allocation error 2 in tensor()");
        m[i] -= jlow;
    }
    /* allocate rows */
    for (i=ilow;i<=ihigh;i++)
        for (j=jlow;j<=jhigh;j++)
        {
            m[i][j]=(double *)malloc((size_t) (khigh-klow+1)*sizeof(double));
            if(!m[i][j]) nrerror("allocation error 3 in tensor()");
            m[i][j] -= klow;
        }
    return m;
}

void free_tensor(double ***m,int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh)
{
    int i,j;
    
    for(i=ihigh;i>=ilow;i--) for(j=jhigh;j>=jlow;j--) free((void*) (m[i][j]+klow));
    for(i=ihigh;i>=ilow;i--) free((void*) (m[i]+jlow));
    free((void*) (m+ilow));
}



/*===============*/

double ****fensor(int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh,int llow,int lhigh)
{
    int i,j,k;
    double ****m;
    
    /* Allocate pointers to the tensors  */
    m=(double ****)malloc((size_t) (ihigh-ilow+1)*sizeof(double***));
    if (!m) nrerror("allocation error 1 in fensor()");
    m -= ilow;
    
    /* allocate pointers to the matrices */
    for (i=ilow;i<=ihigh;i++)
    {
        m[i]=(double ***)malloc((size_t) (jhigh-jlow+1)*sizeof(double**));
        if(!m[i]) nrerror("allocation error 2 in fensor()");
        m[i] -= jlow;
    }
    /* allocate pointers to the rows */
    for (i=ilow;i<=ihigh;i++)
        for (j=jlow;j<=jhigh;j++)
        {
            m[i][j]=(double **)malloc((size_t) (khigh-klow+1)*sizeof(double*));
            if(!m[i][j]) nrerror("allocation error 3 in fensor()");
            m[i][j] -= klow;
        }
    /* allocate rows */
    for (i=ilow;i<=ihigh;i++)
        for (j=jlow;j<=jhigh;j++)
            for(k=klow;k<=khigh;k++)
            {
                m[i][j][k]=(double *)malloc((size_t) (lhigh-llow+1)*sizeof(double));
                if(!m[i][j][k]) nrerror("allocation error 4 in fensor()");
                m[i][j][k] -= llow;
            }
    return m;
}

void free_fensor(double ****m,int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh,int llow,int lhigh)
{
    int i,j,k;
    
    for(i=ihigh;i>=ilow;i--) for(j=jhigh;j>=jlow;j--) for(k=khigh;k>=klow;k--) free((void*) (m[i][j][k]+llow));
    for(i=ihigh;i>=ilow;i--) for(j=jhigh;j>=jlow;j--) free((void*) (m[i][j]+klow));
    for(i=ihigh;i>=ilow;i--) free((void*) (m[i]+jlow));
    free((void*) (m+ilow));
}

/*===============*/





int ***itensor(int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh)
{
    int i,j;
    int ***m;
    
    /* Allocate pointers to the matrices  */
    m=(int ***)malloc((size_t) (ihigh-ilow+1)*sizeof(int**));
    if (!m) nrerror("allocation error 1 in itensor()");
    m -= ilow;
    
    /* allocate pointers to the rows */
    for (i=ilow;i<=ihigh;i++)
    {
        m[i]=(int **)malloc((size_t) (jhigh-jlow+1)*sizeof(int*));
        if(!m[i]) nrerror("allocation error 2 in itensor()");
        m[i] -= jlow;
    }
    for (i=ilow;i<=ihigh;i++)
        for (j=jlow;j<=jhigh;j++)
        {
            m[i][j]=(int *)malloc((size_t) (khigh-klow+1)*sizeof(int));
            if(!m[i][j]) nrerror("allocation error 3 in itensor()");
            m[i][j] -= klow;
        }
    return m;
}

void free_itensor(int ***m,int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh)
{
    int i,j;
    
    for(i=ihigh;i>=ilow;i--) for(j=jhigh;j>=jlow;j--) free((void*) (m[i][j]+klow));
    for(i=ihigh;i>=ilow;i--) free((void*) (m[i]+jlow));
    free((void*) (m+ilow));
}










//================================
//================================



double erfcc(double x)
//Returns the complementary error function erfc(x)
//with fractional error everywhere less than 1.2e-7
{
    double t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
                                                             t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                                                                                                            t*(-0.82215223+t*0.17087277)))))))));
    return x >= 0.0 ? ans : 2.0-ans;
}


//=============================

static int wi[9]={3956, 23552, -3712, 41984, -18160, 41984, -3712, 23552, 3956};
static int LEVMAX = 10;

double quad8stp(double (* fun)(double,double *),double a,double b,double tol,int lev,double *w,double *x0,double *f0,double Q0,int *precur_lev_excess,double *OPT){
    int i;
    double x[17],f[17],Q1,Q2,Q,h,c;
    
    for(i=0;i<17;i++) x[i]=f[i]=0.0;
    for(i=0;i<17;i+=2) x[i]=x0[i/2],f[i]=f0[i/2];
    for(i=0;i<16;i+=2) x[i+1]=(x0[i/2]+x0[i/2+1])/2,f[i+1]=(*fun)(x[i+1],OPT);
    
    h = (b-a)/16;
    for(i=0,Q1=Q2=0.0;i<9;i++) Q1+=h*w[i]*f[i],Q2+=h*w[i]*f[i+8];
    Q = Q1 + Q2;
    
    if (fabs(Q - Q0) > tol*fabs(Q)  && lev <= LEVMAX) {
        c = (a+b)/2;
        Q1= quad8stp(fun,a,c,tol/2,lev+1,w,x,f,Q1,precur_lev_excess,OPT);
        Q2= quad8stp(fun,c,b,tol/2,lev+1,w,x+8,f+8,Q2,precur_lev_excess,OPT);
        Q = Q1 + Q2;
    }
    else if (lev > LEVMAX)
        if (!(*precur_lev_excess)) {
            printf("Recursion level limit reached in quad8. Singularity likely.\n");
            *precur_lev_excess = 1;
        }
        else
            (*precur_lev_excess) ++;
    
    return Q;
}



double quad8(double (* fun)(double,double *),double a,double b,double tol,double *OPT) {
    int i,recur_lev_excess=0;
    double x[9],y[9],w[9],Q;
    
    // Top level initialization, Newton-Cotes weights
    for(i=0;i<9;i++) {
        w[i] = wi[i]/14175.0;
        x[i] = a + i*(b-a)/8.0;
        y[i] = (*fun)(x[i],OPT);
    }
    //Adaptive, recursive Newton-Cotes 8 panel quadrature
    Q= quad8stp(fun,a,b,tol,1,w,x,y,1.0e100,&recur_lev_excess,OPT);
    if (recur_lev_excess > 1) printf("Recursion level limit reached %d times.\n",recur_lev_excess);
    
    return Q;
}




void quicksort(double *p, int left, int right) {
    int i, j;
    double pivot;
    
    if (left < right) {
        i = left;
        j = right + 1;
        pivot = p[left];
        do {
            do
                i++;
            while ((p[i] < pivot) && (i < right));
            do
                j--;
            while ((p[j] > pivot) && (j > left));
            if (i < j) {
                SWAP(&p[i],&p[j]);
            }
        } while (i < j);
        SWAP(&p[left], &p[j]);
        quicksort(p, left, j-1);
        quicksort(p, j+1, right);
    }
    
}


void quicksort1(double *p, int *index, int left, int right) {
    int i, j;
    double pivot;
    
    
    if (left < right) {
        i = left;
        j = right + 1;
        pivot = p[left];
        do {
            do 
                i++;
            while ((p[i] < pivot) && (i < right));
            do 
                j--;
            while ((p[j] > pivot) && (j > left));
            if (i < j) {
                SWAP(&p[i],&p[j]);
                SWAPi(&index[i],&index[j]);
            }
        } while (i < j);
        SWAP(&p[left], &p[j]);
        SWAPi(&index[left], &index[j]);
        quicksort1(p,index, left, j-1);
        quicksort1(p,index, j+1, right);
    }
    
}



void SWAP(double *a, double *b) {
    double temp;
    
    temp=*a;
    *a=*b;
    *b=temp;
}

void SWAPi(int *a, int *b) {
    int temp;
    
    temp=*a;
    *a=*b;
    *b=temp;
}