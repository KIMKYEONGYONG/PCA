// CYPCA.cpp: implementation of the CCYPCA class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "CYPCA.h"
#include <math.h>
#include <iostream>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// 아래에서 사용된 HH 표기 = HouseHolder의 HH



// 하우스홀더(HouseHolder) Method
void CCYPCA::HH_Eigen(CCYMatrix &vector, CCYMatrix &value, CCYMatrix &A)
{
	value.Create(1, A.Col());
	vector.Create(A.Row(), A.Col());

	double **s = new double *[A.Row() + 1];   ///////////////////
	int r;
	for (r = 0; r <= A.Row(); r++) {
		s[r] = new double[A.Col() + 1];       ///////////////////
	}
	double *sp = new double[A.Row() + 1];     ///////////////////
	double *e = new double[A.Row() + 1];      ///////////////////

	int n, i, j;
	n = A.Row();

	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++) {
			s[i + 1][j + 1] = A[i][j];
		}
	}

	_HH_tred2(s, n, sp, e);
	_HH_tqli(sp, e, n, s);

	for (i = 0; i<n; i++) {
		for (j = 0; j<n; j++) {
			vector[i][j] = s[i + 1][j + 1];
		}
		value[0][i] = sp[i + 1];
	}

	for (r = 0; r <= A.Row(); r++) {
		delete[] s[r];
	}

	delete [] e;
	delete [] s;
	delete [] sp;
}


double CCYPCA::_HH_pythag(double a, double b)
{
	double absa, absb;
	absa = fabs(a);
	absb = fabs(b);

	if (absa > absb) return absa * sqrt(1.0 + POWER(absb / absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + POWER(absa / absb)));
}


void CCYPCA::_HH_tqli(double d[], double e[], int n, double **z)
{
	int m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;
	for (i = 2; i <= n; i++) e[i - 1] = e[i]; // Convenient to renumber the elements of e.
	e[n] = 0.0;
	
	for (l = 1; l <= n; l++)
	{
		iter = 0;
		do  // m != l
		{
			for (m = l; m <= n - 1; m++) {
				//Look for a single small subdiagonal element to split the matrix.
				dd = fabs(d[m]) + fabs(d[m + 1]);
				//if ((float)(fabs(e[m]) + dd) == dd) break;
				if ((fabs(e[m]) + dd) == dd) break;
			}
			

			if (m != l){

				//if (iter++ == 30) { std::cout << "\ttqli() : Too many iterations in tqli"; }

				iter++;

				g = (d[l + 1] - d[l]) / (2.0*e[l]); //Form shift.
				r = _HH_pythag(g, 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g)); //This is dm . ks.
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--)
				{
					//QL, followed by Givens rotations to restore tridiagonal form.
					f = s*e[i];
					b = c*e[i];
					e[i + 1] = (r = _HH_pythag(f, g));
					if (r == 0.0)
					{
						//Recover from underflow.
						d[i + 1] -= p;
						e[m] = 0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g)*s + 2.0*c*b;
					d[i + 1] = g + (p = s*r);
					g = c*r - b;
					/* Next loop can be omitted if eigenvectors not wanted*/
					for (k = 1; k <= n; k++)
					{
						//Form eigenvectors.
						f = z[k][i + 1];
						z[k][i + 1] = s*z[k][i] + c*f;
						z[k][i] = c*z[k][i] - s*f;
					}

				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
			} // if (m != l )

		} while (m != l);
	}



}

void CCYPCA::_HH_tred2(double **a, int n, double d[], double e[])
{
	int l, k, j, i;
	double scale, hh, h, g, f;
	for (i = n; i >= 2; i--)
	{
		l = i - 1;
		h = scale = 0.0;
		if (l > 1)
		{
			for (k = 1; k <= l; k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0) //Skip transformation.
				e[i] = a[i][l];
			else
			{
				for (k = 1; k <= l; k++)
				{
					a[i][k] /= scale; //Use scaled a for transformation.
					h += a[i][k] * a[i][k]; //Form s in h.
				}
				f = a[i][l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale*g;
				h -= f*g; //Now h is equation (11.2.4).
				a[i][l] = f - g; //Store u in the ith row of a.
				f = 0.0;
				for (j = 1; j <= l; j++)
				{
					/* Next statement can be omitted if eigenvectors not wanted */
					a[j][i] = a[i][j] / h; //Store u/H in ith column of a.
					g = 0.0; //Form an element of A ?u in g.
					for (k = 1; k <= j; k++)
						g += a[j][k] * a[i][k];
					for (k = j + 1; k <= l; k++)
						g += a[k][j] * a[i][k];
					e[j] = g / h; //Form element of p in temporarily unused element of e.
					f += e[j] * a[i][j];
				}
				hh = f / (h + h);
				for (j = 1; j <= l; j++)
				{ //Form q and store in e overwriting p.
					f = a[i][j];
					e[j] = g = e[j] - hh*f;
					for (k = 1; k <= j; k++) //Reduce a, equation (11.2.13).
						a[j][k] -= (f*e[k] + g*a[i][k]);
				}
			}
		}
		else
			e[i] = a[i][l];
		d[i] = h;
	}
	/* Next statement can be omitted if eigenvectors not wanted */
	d[1] = 0.0;
	e[1] = 0.0;
	/* Contents of this loop can be omitted if eigenvectors not
	wanted except for statement d[i]=a[i][i]; */
	for (i = 1; i <= n; i++)
	{ //Begin accumulation of transformation matrices.
		l = i - 1;
		if (d[i]) { //This block skipped when i=1.
			for (j = 1; j <= l; j++)
			{
				g = 0.0;
				for (k = 1; k <= l; k++) //Use u and u/H stored in a to form P·Q.
					g += a[i][k] * a[k][j];
				for (k = 1; k <= l; k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i] = a[i][i]; //This statement remains.
		a[i][i] = 1.0; //Reset row and column of a to identity matrix for next iteration. 
		for (j = 1; j <= l; j++) a[j][i] = a[i][j] = 0.0;
	}

}