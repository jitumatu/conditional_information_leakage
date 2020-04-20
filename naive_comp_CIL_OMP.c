/* To compile:  $ gcc -fopenmp -O3 <filename>.c -lm -o <filename>.exe */

#include <stdio.h>
#include <stdlib.h>   /* for atoi */
#include <omp.h>      /* for parallel processing */
#include <time.h>     /* for clock_gettime */
#include <math.h>
#include <limits.h>

#define X 2 // size of input alphabet
#define Z 2 // size of output alphabet 

#define n 30 // code length
#define m 10 // message size is 2^m 

#define delta 0.1 /* cross-over probability */


int main(){
	struct timespec startTime, endTime; /* count the time */ 
	double W[Z][X]; /* transition probability matrix  of the channel */
	double Phi[X][Z]; /* transition probability matrix of backward chanel */
	double Phi0, Phi1;
	double H=0, temp, *beta, *beta_new, Pr;
	unsigned long long int a[n]; /* integer representation of n-th column of parity check matrix A */ 
	int z[n]; // the output symbol.
	int MAXthread;
	long long int i, rnd;
	int j, k, p, r, s, t, v, x, y, num;
	
	MAXthread = omp_get_max_threads(); 
	printf("max threads (set): %d\n", MAXthread);

	/* we have to generate a random number larger than RAND_MAX*/
	printf("RAND_MAX=%d", RAND_MAX);
	num = RAND_MAX;
	k=0;
	while (num !=0){
		num = num/2;
		k++;
	}
	p = m/k; /* p+1 is the number of rand() to generate random number up to 2^m-1 */
	printf("p=%d\n", p);
	
	/* Transition probability of Eve's channel W_E(z|x) */
	/* Not necessarily a BSC */
	W[0][0] = W[1][1] = 1-delta;
	W[0][1] = W[1][0] = delta;
		
	/* Determine backward channel. 
	   We use y to denote the symbol z, 
	   because z is used for the received signal */
	for (y = 0; y < Z; y++){
		temp = 0;
		for (x = 0; x < X; x++){
			temp += W[y][x];
		}
		for (x = 0; x < X; x++){
			Phi[x][y] = W[y][x] /= temp;
			printf("# Phi(%d|%d) = %f\n", x, y, Phi[x][y] );
		}		
	}
	
	/* determine the parity matrix A */
	/* 1<< m is "2 to the m-th power" */
	
	for (i = 0; i < m; i++){
		a[i] = 1<<i;
		printf("a[%lld]=%ld\n", i, a[i]);
	}
	
	for (i = m; i < n; i++){
		rnd = rand();
		printf("rnd=%ld\n", rnd);
		for (j=0; j<p; p++);
			rnd = rnd * RAND_MAX + rand();
			a[i] = rnd % ( 1<<m );
		printf("a[%lld]=%ld\n", i, a[i]);
	}
	
	/* deterine the reseived symbols */
	for (i = 0; i < n; i++){
		z[i] = rand()%2;
		// printf("z[%d]=%d\n", i, z[i]);
	}

	
	beta = (double *) malloc(sizeof(double) * ( ( (long long int) 1)  << m)); 

	for( i = 0; i < ( ( (long long int) 1) << m); i++ ){
		beta[i] = 0.0;
	}
	

	/* get start_time */
	clock_gettime(CLOCK_REALTIME, &startTime);

	#pragma omp parallel for private(s, t, r, v, Pr) reduction(+:beta[0:1<<m]) 
	for (i=0; i< ( ( (long long int) 1) << n ); i++){
		s = 0; 
		t = 0;
		Pr = 1;
		for (r = 0; r < n; r++){
			v = ( ( i & ( (long long int) 1) << r ) == 0 ) ? 0 : 1;
			s = s ^ ( a[r] * v );
			Pr = Pr * Phi[ v ][ z[r] ];
		}
		beta[s] += Pr;
	}
	
	/* get end_time */
	clock_gettime(CLOCK_REALTIME, &endTime);

	for (i = 0; i < ( ( (long long int) 1) << m ) ; i++){
		if (beta[i] != 0){
			H += - beta[i] * log( beta[i] );
		}
		printf("beta[%lld]=%f\n", i, beta[i]);
	}

	H = H/log(2); /* the base of log is 2 */
	printf("I(S^m;Z^n) = %f\n", m - H);
	
	
	/* print the elapsed time */
	printf("elapsed time = ");
	if (endTime.tv_nsec < startTime.tv_nsec) {
		printf("%5ld.%09ld", endTime.tv_sec - startTime.tv_sec - 1,
		       endTime.tv_nsec + (long int)1.0e+9 - startTime.tv_nsec);
	} else {
		printf("%5ld.%09ld", endTime.tv_sec - startTime.tv_sec,
		       endTime.tv_nsec - startTime.tv_nsec);
	}
	printf("(sec)\n");
	
	free(beta);
	free(beta_new);
	return 0;
}
