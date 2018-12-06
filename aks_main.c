/**
* @file aks_main.c
* the main file of the aks-project (GMP-Version)
* @author Alexandru Paler and Fabio Campos
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "gmp.h"
#include "newton_gmp.h"
#include "findR.h"
#include "euler.h"
#include "polyUtil.h"

/* timing-var. */
struct timeval timev1, timev2; 

/**
 * Main entry point of the program
 * @param	number of arguments
 * @return 	0 if everyting ok
*/
int main(int argc, char* argv[]) {

	int i;
	int logn;
	int myLimit;
	mpz_t a;
	mpz_t myR;
	int myRint, myEulerCounter;
	
	mpz_t dummy;
	int isPrime;
	int myBinomLimit;
	
	poly basisPoly;
	poly rightPoly;
	poly resultPoly;

	mpz_init(a);
	mpz_init(myR);

	if (mpz_set_str(a, argv[1], 0) != 0) {
		/* this is a prototype without error handling */
	}

	gmp_printf("%Zd,", a);
	gmp_printf("1,");

	gettimeofday(&timev1, NULL);

	if ((newton_it(a))) {
		(void) printf("composite\n");
		return 0;
		/* it is composite */
	}
	gettimeofday(&timev2, NULL);

	(void) printf("%f,", (timev2.tv_sec - timev1.tv_sec) + 0.000001
			* (timev2.tv_usec - timev1.tv_usec));


	gettimeofday(&timev1, NULL);

	myRint=findR(a);

	(void) printf("%d,", myRint);
	if (myRint==-1) {
		(void) printf("prime\n");
		return 0;
		/* it is composite */
	}

	gettimeofday(&timev2, NULL);

	(void) printf("%f\n", (timev2.tv_sec - timev1.tv_sec) + 0.000001
			* (timev2.tv_usec - timev1.tv_usec));
	//return 0;
	gettimeofday(&timev1, NULL);

	myEulerCounter = euler_probDiv(myRint);

	logn = (int)mpz_sizeinbase(a, 10);

	myLimit = sqrt(myEulerCounter) * logn;

	mpz_init(dummy);

	basisPoly.degree = 2;

	mpz_mod_ui(dummy, a, myRint);

	rightPoly.degree = mpz_get_ui(dummy) + 1;

	/* basisPoly => (x + a) */
	basisPoly.coeff = (mpz_t*)malloc(basisPoly.degree * sizeof(mpz_t));
	polyInit(&basisPoly);
	mpz_set_ui(basisPoly.coeff[1], 1);

	/* rightPoly => (x^(n mod r) + a) */
	rightPoly.coeff = (mpz_t*)malloc(rightPoly.degree * sizeof(mpz_t));
	polyInit(&rightPoly);
	mpz_set_ui(rightPoly.coeff[rightPoly.degree - 1], 1);

	/* the bit-length of the input */
	myBinomLimit = mpz_sizeinbase(a, 2);

	isPrime=1;

	gettimeofday(&timev1, NULL);

	fflush(stdout);

	for (i=1; i<= myLimit; i++) {

		mpz_set_ui(basisPoly.coeff[0], i);
		mpz_set_ui(rightPoly.coeff[0], i);

		/* calculate (x + a)^n mod (x^r - 1, n) */
		polyExp_Mod_int(&resultPoly, &basisPoly, a, myRint, myBinomLimit);

		/*
		 * check if (x + a)^n mod (x^r - 1) != x^(n mod r) + a
		 */
		if (polyIsEqual(&resultPoly, &rightPoly)!=1) {
			printf("composite,");
			isPrime=0;
			break;
		}
	}

	gettimeofday(&timev2, NULL);

	if (isPrime==1) {
		printf("prime,");
	}

	(void) printf("%f\n", (timev2.tv_sec - timev1.tv_sec) + 0.000001
			* (timev2.tv_usec - timev1.tv_usec));

	return 0; /* everything is allright */

}
