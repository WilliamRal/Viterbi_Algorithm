//
//  utils.h
//  Project1_555
//
//  Created by Bing Hong Fu on 2/2/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//

#ifndef __Project1_555__utils__
#define __Project1_555__utils__

#include <stdio.h>
int uni_int(int *idum,int m);

int uni_int2(long int *idum,int m);

double gauss1(int *idum);

double gauss2(long int *idum);

double ran1(int *idum);

double ran2(long int *idum);

void nrerror(char *errortext);

int ipow(int m, int p);

double mag_sq(double *v,int ilow,int ihigh);



int *ivector(int ilow,int ihigh);

void free_ivector(int *v,int ilow,int ihigh);

int **imatrix(int ilow,int ihigh,int jlow,int jhigh);

void free_imatrix(int **m,int ilow,int ihigh,int jlow,int jhigh);

char **charmatrix(int ilow,int ihigh,int jlow,int jhigh);

void free_charmatrix(char **m,int ilow,int ihigh,int jlow,int jhigh);

int ***itensor(int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh);

void free_itensor(int ***m,int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh);



double *vector(int ilow,int ihigh);

void free_vector(double *v,int ilow,int ihigh);

double **matrix(int ilow,int ihigh,int jlow,int jhigh);

void free_matrix(double **m,int ilow,int ihigh,int jlow,int jhigh);

double ***tensor(int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh);

void free_tensor(double ***m,int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh);

double ****fensor(int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh,int llow,int lhigh);

void free_fensor(double ****m,int ilow,int ihigh,int jlow,int jhigh,int klow,int khigh,int llow,int lhigh);



double erfcc(double x);

double quad8(double (* fun)(double,double *),double a,double b,double tol,double *OPT);

void quicksort(double *p, int left, int right);

void quicksort1(double *p,int *index, int left, int right);

void SWAP(double *a, double *b);

void SWAPi(int *a, int *b);

#endif /* defined(__Project1_555__utils__) */
