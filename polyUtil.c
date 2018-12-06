/**
 * @file polyUtil.c
 * functions to manipulate the polynom-struct
 * @author Alexandru Paler and Fabio Campos
 * @project AKS on Cell - WS2007_08@FH-Wiesbaden
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "polyUtil.h"

/**
 * Prints out an instance of the poly-struct
 * @param  	   	inPoly	a reference to the poly that should be printed
 * @return      void
 */
void polyPrint(poly* inPoly) {
	int i;
	(void) printf("(");
	for (i=inPoly->degree-1; i>=0; i--) {
		(void) gmp_printf(" %Zdx^%d ", inPoly->coeff[i], i);
	}
	(void) printf(")\n");
}

/**
 * Initializes the elements of a poly-struct
 * @param  	   	inPoly	a reference to the poly-struct
 * @return      void
 */
void polyInit(poly* inPoly) {
	int i;

	for (i=0; i<= (inPoly->degree-1); i++) {
		mpz_init(inPoly->coeff[i]);
	}
}

/**
 * Makes a deep copy of the input struct
 * @param  	   	resultPoly	a reference to the result polynom 
 * @param  	   	inPoly		a reference to the polynom, which shoud be copied
 * @return      the result polynom
 */
inline void polyCopy(poly* resultPoly, poly* inPoly) {
	int i;

	resultPoly->degree = inPoly->degree;

	resultPoly->coeff = (mpz_t *)malloc(resultPoly->degree * sizeof(mpz_t));

	for (i=0; i<= (inPoly->degree-1); i++) {

		mpz_init(resultPoly->coeff[i]);

		mpz_set(resultPoly->coeff[i], inPoly->coeff[i]);

	}

}

/**
 * Compares two poly-structs
 * @param  	   	left	a reference to the first polynom
 * @param  	   	right	a reference to the second polynom
 * @return      the result polynom
 */
int polyIsEqual(poly* left, poly* right) {
	int i;

	for (i=0; i<=left->degree; i++) {
		if ((mpz_cmp_ui(left->coeff[i], 0)!=0) && (mpz_cmp(left->coeff[i],
				right->coeff[i])!=0)) {
			return 0;
		}
	}
	return 1;
}

/**
 * This function calculates the exponentiation of a polynom-struct using
 * the repeated squaring method
 * @param  	   	resultPoly	a reference to the result polynom
 * @param		left		a reference to the polynom that should be calculated
 * @param		exp			the exponent
 * @param		r			the r-value
 * @param		myLimit		the bit-length of exp
 * @return      the result polynom
 */
void polyExp_Mod_int(poly* resultPoly, poly* left, mpz_t exp, int r, int myLimit) {
	poly* dummyPoly;

	int k;

	dummyPoly = (poly*) malloc(sizeof(poly));

	polyCopy(resultPoly, left);

	for (k=myLimit-2; k>-1; k--) {

		polyMult_Mod(dummyPoly, resultPoly, resultPoly, exp, r);

		/* checks if the i-th bit == 1 */
		if (mpz_tstbit(exp, k)==1) {

			polyMult_Mod(resultPoly, dummyPoly, left, exp, r);

		} else {

			polyCopy(resultPoly, dummyPoly);

		}

	}
	free(dummyPoly->coeff);
}

/**
 * This function calculates the multiplication of two polynom-structs
 * and cut the result-polynom by mod (x^r - 1, exp)
 * @param  	   	resultPoly	a reference to the result polynom
 * @param		left		a reference to the polynom that should be multiplied
 * @param		right		a reference to the polynom that should be multiplied
 * @param		exp			the mod-value
 * @param		r-value		the r-value (x^r - 1)
 * @return      the result polynom
 */
inline void polyMult_Mod(poly* resultPoly, poly* left, poly* right, mpz_t exp,
		int r) {
	int i, j;
	int countMod;
	mpz_t dummy;

	mpz_init(dummy);

	if ((left->degree + right->degree -1) < r) {
		resultPoly->degree = left->degree + right->degree -1;
	} else {
		resultPoly->degree = r;
	}

	resultPoly->coeff = (mpz_t *)malloc(resultPoly->degree * sizeof(mpz_t));

	for (i=0; i<= (resultPoly->degree-1); i++) {
		mpz_init(resultPoly->coeff[i]);
	}

	for (i=0; i<= (left->degree-1); i++) {
		for (j=0; j<= (right->degree-1); j++) {
			countMod = (i+j)%r;
			mpz_mul(dummy, left->coeff[i], right->coeff[j]);
			mpz_add(resultPoly->coeff[countMod], resultPoly->coeff[countMod],
					dummy);
			mpz_mod(resultPoly->coeff[countMod], resultPoly->coeff[countMod],
					exp);
		}
	}

}
