/*
  #############################################################
  ##   BSPedupack Version 2.0                                ##
  ##   Copyright (C) 2016 Rob H. Bisseling                   ##
  ##                                                         ##
  ##   BSPedupack is released under the                      ##
  ##   GNU LESSER GENERAL PUBLIC LICENSE                     ##
  ##   Version 3, 29 June 2007 (given in the file LICENSE)   ##
  #############################################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <complex.h>
#include <math.h>
#include <tgmath.h>

#include "mcbsp.h" // header file for MulticoreBSP for C

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

typedef struct item {
  double weight ;
  long index ;
} Item ;

long *vecalloci(size_t n);
double *vecallocd(size_t n);
complex *vecallocc(size_t n);
Item *vecallocitem(size_t n);
double **matallocd(size_t m, size_t n);
void vecfreei(long *pi);
void vecfreed(double *pd);
void vecfreec(complex *pc);
void vecfreeitem(Item *pitem);
void matfreed(double **ppd);
