#include "bspedupack.h"
#include "sys/time.h"
#include <stdlib.h>
#include <time.h>

long P;
long N;
long x;
long y;

void matmul(){
	bsp_begin(P);
	long p = bsp_nprocs();
	long s = bsp_pid();
	long k = N/p;
	double *A;
	double *B;
	double *C;
	double **pA;
	double **pB;
	double **pC;
	struct timeval t1;
	int i,j,g;

	A = (double *) malloc (N*k*sizeof(double));
	B = (double *) malloc (N*k*sizeof(double));
	C = (double *) malloc (N*k*sizeof(double));
	
	pA = (double **) malloc(k*sizeof(double *));
	pB = (double **) malloc(k*sizeof(double *));
	pC = (double **) malloc(k*sizeof(double *));

	for(i=0; i<k; i++){
		pA[i]=A+i*N;
		pB[i]=B+i*N;
		pC[i]=C+i*N;
	}

	gettimeofday(&t1, NULL);

	srand(t1.tv_usec * t1.tv_sec * s);
	for(i=0; i<k; i++)
		for(j=0; j<N; j++)
			pA[i][j] = (double)rand()/(double)(RAND_MAX);

	for(i=0; i<k; i++)
		for(j=0; j<N; j++){
			pB[i][j] = (double)rand()/(double)(RAND_MAX);
			pC[i][j] = 0.0;
		}

/* ----- Main algorithm -----*/
	
   	bsp_push_reg(B, N*k*sizeof(double));
	bsp_sync();

	
	int l;
	int startingpoint = s*k;
	int endpoint = startingpoint + k;
	int row;
	int superstep;
	double tp1, tp2=0;
	double tp3, tp4=0;
	double time0= bsp_time();
	
	for(superstep = 0; superstep<p; superstep++){
	
		tp3 = bsp_time();
		for (row=0; row<k ; row++)
			for (j=startingpoint, l=0; j<endpoint; j++, l++)
				for (g=0; g<N; g++)
					pC[row][g] += pA[row][j] * pB[l][g];				
		
		tp4 += bsp_time()-tp3;

		tp1 = bsp_time();
		bsp_get((s+1)%p, B, 0, B, N*k*sizeof(double) );
		bsp_sync();
		tp2 += bsp_time()-tp1;
		
		startingpoint = (startingpoint+k)%N;
		endpoint = startingpoint +k ;
	}
	
	double time1= bsp_time();
/* ----- Print main algorithm result-----*/
	
	if (s==y/k){
		printf("---- C[%ld][%ld]=%.2f from row-wise decomposition-----\n", y, x, pC[y%k][x]);
		printf("Total time : %f, communication time: %f, computation time %f: \n", time1-time0, tp2, tp4);
	}	
/* ----- test -----*/

	double *bcolumn;
	bcolumn = (double *) malloc (N*sizeof(double));
	bsp_push_reg(bcolumn, N*sizeof(double));
	bsp_sync();

	for(i=0; i<k; i++){
		bsp_put(y/k, &(pB[i][x]), bcolumn, (s*k+i)*sizeof(double) , sizeof(double) );			
		bsp_sync();
	}
	if(s==y/k){		
		double sum=0.0;
		for(i=0; i<N; i++)
			sum += pA[y%k][i]*bcolumn[i];	
		printf("---- C[%ld][%ld]=%.2f from test -----\n", y, x, sum);
		
	}	

	bsp_pop_reg(bcolumn);
	bsp_pop_reg(B);

	free(A);
	free(B);
	free(C);
	free(bcolumn);
	free(pA);
	free(pB);
	free(pC);
	bsp_sync();
	bsp_end();
	return;
}

int main(int argc, char **argv){
	bsp_init(matmul, argc, argv);

    printf("how many processors do you want to use ?\n");
    scanf("%ld",&P);
    if (P>bsp_nprocs()){
        printf("sorry, only %u processors are available\n", bsp_nprocs());
		return -1;
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

