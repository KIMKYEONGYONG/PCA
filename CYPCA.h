// CYPCA.h: interface for the CCYPCA class.
//
//////////////////////////////////////////////////////////////////////

#ifndef __CCYPCA_CLASS_HEADER__
#define __CCYPCA_CLASS_HEADER__

/*
PCA(Principal Component Analysis)
1) 하우스홀더 알고리즘을 사용한 고유치, 고유벡터 산출

2005. 12. by Chi-Young Cho.
*/

#include "stdafx.h"
#include "CYMatrix.h"
#define POWER(a) ((a)*(a))
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

class CCYPCA
{
private:
	/****************************************************/
	/* Householder 방법을 이용한 고유치/고유벡터 산출   */
	/****************************************************/
	static double _HH_pythag(double a, double b);
	static void _HH_tred2(double **a, int n, double d[], double e[]);
	static void _HH_tqli(double d[], double e[], int n, double **z);

public:
	static void HH_Eigen(CCYMatrix &vector, CCYMatrix &value, CCYMatrix &A); // householder 메소드를 이용

};

#endif