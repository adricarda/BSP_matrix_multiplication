#include "bspedupack.h"
#include "sys/time.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

long P;
long N;
long x;
long y;

void localMatUpdate(double **pA, double **pB, double **pC, long dim){
	int i,j,k;
    for (i=0;i<dim;i++) {
      for (k=0;k<dim;k++) {
        for (j=0;j<dim;j++) {
          pC[i][j] += pA[i][k]*pB[k][j];
        }
      }
    }
}

void matmul(){
	bsp_begin(P);
	long p = bsp_nprocs();
	long s = bsp_pid();
	long pside = sqrt(p);
	long k = N/pside;
	long row;
	long column;
	double *A, *B, *C;
	double **pA, **pB, **pC;
	struct timeval t1;
	int i,j;

	A = (double *) malloc (k*k*sizeof(double));
	B = (double *) malloc (k*k*sizeof(double));
	C = (double *) malloc (k*k*sizeof(double));
	
	pA = (double **) malloc(k*sizeof(double *));
	pB = (double **) malloc(k*sizeof(double *));
	pC = (double **) malloc(k*sizeof(double *));

	row = s/pside;
	column = s%pside;

	long prow = y/k;
	long pcol = x/k;
	long id = prow*pside+pcol;

	bsp_push_reg(A, k*k*sizeof(double));
	bsp_push_reg(B, k*k*sizeof(double));
	bsp_push_reg(C, k*k*sizeof(double));	
	bsp_sync();	

	for(i=0; i<k; i++){
		pA[i]=A+i*k;
		pB[i]=B+i*k;
		pC[i]=C+i*k;
	}

	gettimeofday(&t1, NULL);

	srand(t1.tv_usec * t1.tv_sec * s);
	for(i=0; i<k; i++)
		for(j=0; j<k; j++)
			pA[i][j] = (double)rand()/(double)(RAND_MAX);

	for(i=0; i<k; i++)
		for(j=0; j<k; j++){
			pB[i][j] = (double)rand()/(double)(RAND_MAX);
			pC[i][j] = 0.0;
		}

/*----------- compute C[y][x] explicity ---------*/ 
	double *arow;
	double *bcolumn;
	arow = (double *) malloc(N*sizeof(double));
	bcolumn = (double *) malloc(N*sizeof(double *));
	bsp_push_reg(arow,N*sizeof(double *));
	bsp_push_reg(bcolumn,N*sizeof(double *));
	bsp_sync();

	if(row == prow){	
			bsp_put(0, pA[y%k] , arow, k*column*sizeof(double), k*sizeof(double));	
	}
	bsp_sync();

	if(column == pcol){	
		for(i=0; i<k; i++){
			bsp_put(0, &pB[i][x%k] , bcolumn, (k*row+i)*sizeof(double) , sizeof(double));	
		}
	}
	bsp_sync();

	if(s==0){
	     double sum=0.0;
        for(i=0; i<N; i++)
            sum += arow[i]*bcolumn[i];
        printf("---- C[%ld][%ld]=%.2f from test -----\n", y, x, sum);
		
	}
   	bsp_pop_reg(arow);
	bsp_pop_reg(bcolumn);
	free(arow);
	free(bcolumn);

/*----------- shifting matrix --------- */

	for (i=0; i< pside; i++){
		if (row > 0 && i<row){
			if (s%pside == pside-1){
				//i'm the last on my row
				bsp_get(s-pside+1, A, 0, A, k*k*sizeof(double));	
			}
			else{
				bsp_get(s+1, A, 0, A, k*k*sizeof(double));	
			}
		}
		bsp_sync();
	}

	for (i=0; i< pside; i++){
		if ( column>0 && i<column){
			if (row == pside-1){
				//i'm the last on my column
				bsp_get((s+pside)%pside, B, 0, B, k*k*sizeof(double));	
			}
			else{
				bsp_get(s+pside, B, 0, B, k*k*sizeof(double));	
			}
		}
		bsp_sync();
	}

/*----------- compute matrix ---------*/ 
	
	double start = bsp_time();	
	for(i=0; i<pside; i++){
		localMatUpdate(pA, pB, pC, k);

		//get matrix A from right
		if (s%pside == pside-1){
			//i'm the last on my row
			bsp_get(s-pside+1, A, 0, A, k*k*sizeof(double));	
		}
		else{
			bsp_get(s+1, A, 0, A, k*k*sizeof(double));	
		}
		//get matrix B from down
		if (row == pside-1){
				//i'm the last on my column
				bsp_get((s+pside)%pside, B, 0, B, k*k*sizeof(double));	
			}
			else{
				bsp_get(s+pside, B, 0, B, k*k*sizeof(double));	
			}
		bsp_sync();
	}
	

	double finish = bsp_time();	

/*----------- show requested element ---------*/ 
	if (s==id){
		printf("element C[%ld][%ld]=%.2f in %f\n", y, x, pC[y%k][x%k], finish-start);
	}

	bsp_pop_reg(A);
	bsp_pop_reg(B);
	bsp_pop_reg(C);
	free(A);
	free(B);
	free(C);
	free(pA);
	free(pB);
	free(pC);
	bsp_end();
	return;
}

int main(int argc, char **argv){
	bsp_init(matmul, argc, argv);

    printf("how many processors do you want to use ?\n");
    scanf("%ld",&P);
    if (P>bsp_nprocs()){
        printf("sorry, only %u processors are available\n", bsp_nprocs());
    }
	printf("Insert matrix dimension\n");
	//check if N/P
    scanf("%ld",&N);

	do{
	printf("Insert Y coordinate\n");
    scanf("%ld",&y);
	}while (y< 0 || y>=N);

	do{
	printf("Insert x coordinate\n");
    scanf("%ld",&x);
	}while (x< 0 || x>=N);

	matmul();
    return 1;
}

