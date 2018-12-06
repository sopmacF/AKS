/**
 * @file polyUtil.h
 * functions to manipulate the polynom-struct
 * @author Alexandru Paler and Fabio Campos
 * @project AKS on Cell - WS2007_08@FH-Wiesbaden
 * 
 */

#include "gmp.h"

#ifndef POLYUTIL_H_
#define POLYUTIL_H_

typedef struct  {

	unsigned int degree; /* this value is the real-poly-degree + 1 
	 * for example: the poly x^5 has a degree value of 6
	 * in this struct
	 */
	mpz_t * coeff; /* Array presenting the coefficients of the polynom */

} poly;
 
void polyPrint(poly* in);

void polyInit(poly* inPoly);

void polyExp_Mod_int(poly* resultPoly, poly* left, mpz_t exp, int r, int myLimit);

inline void polyMult_Mod(poly* resultPoly, poly* left, poly* right, mpz_t exp, int r);

int polyIsEqual(poly* left, poly* right);

inline void polyCopy(poly* ResultPoly, poly* inPoly);

int tstbit(int word, int bit);

#endif /*POLYUTIL_H_*/
