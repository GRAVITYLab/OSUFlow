/*
THIS SOFTWARE BELONGS TO THE AMERICAN PEOPLE.  IT WAS DEVELOPED
AT THEIR EXPENSE AT NASA AMES RESEARCH CENTER.

NAME: polynomial solver

REQUIREMENTS: solve various polynomial equations

DISCUSSION:

HISTORY
Originated		????		Creon Levit
*/

/* dont forget that these deal with floats! */

#include "Topology.h"
#include "zero_test.h"

double cube_root ( double x )
{
	if (x > 0.0)
		return (pow (x, 1.0/3.0));
	else if (x < 0.0)
		return (-pow (-x, 1.0/3.0));
	else
		return 0.0;
}

int solve_cubic ( float a, float b, float c, float d,
				 float* x1, float* x2, float* x3 )
{
	/* 
	*         3     2
	* Solve ax  + bx  + cx + d = 0
	*
	* Return values: 0 - no solution 
	*               -1 - infinite number of solutions
	*                1 - only one distinct root, real, possibly of multiplicity 2 or 3.
	*                    returned in *x1
	*		    2 - only two distinct roots, both real, one possibly of multiplicity 2,
	*                    returned in *x1 and *x2
	*		   -2 - only two distinct roots (complex conjugates),
	*		        real part returned in *x1,
	*			imaginary part returned in *x2
	*                3 - three distinct real roots of multiplicity 1
	*                    returned in *x1, *x2, and *x3.
	*               -3 - one real root and one complex conjugate pair.
	*                    real root in *x1, real part of complex conjugate root
	*                    in *x2, imaginary part in *x3.
	*
	*                XXX - this whole scheme is wrong.  
	*                      Each root should be returned along with its multiplicity.
	*/

	double p, q, r, a2, b2, z1, z2, z3, des;

	*x1 = *x2 = *x3 = 0.0;

	if (a == 0.0)
		return (solve_quadratic (b, c, d, x1, x2));
	else
	{	
		p = b/a;			/* reduce to y^3 + py^2 + qy + r = 0 */
		q = c/a;			/* by dividing through by a */
		r = d/a;

		a2 = (3*q - p*p)/3;	        /* reduce to z^3 + a2*z + b2 = 0 */
		Zero_test2 (a2, 3*q, p*p);

		b2 = (2*p*p*p - 9*p*q + 27*r)/27;	/* using y = z - p/3 */
		Zero_test3 (b2, 2*p*p*p, 9*p*q, 27*r);

		des = (b2*b2/4 + a2*a2*a2/27);
		Zero_test2 (des, b2*b2/4, a2*a2*a2/27);

		if (des == 0.0)	/* three real roots, at least two equal */ 
		{
			double a3 = cube_root (-b2/2);
			if (a3 == 0.0)	/* one distinct real root of multiplicity three */
			{
				z1 = 0.0;
				*x1 = z1 - p/3;
				Zero_test2 (*x1, z1, p/3);
				*x2 = 0.0;
				*x3 = 0.0;
				return (1);
			}
			else			/* one real root multiplicity one, another of multiplicity two */
			{
				z1 = 2*a3;
				z2 = -a3;
				*x1 = z1 - p/3;
				Zero_test2 (*x1, z1, p/3);
				*x2 = z2 - p/3;
				Zero_test2 (*x2, z2, p/3);
				*x3 = 0.0;
				return (2);
			}
		}
		else if (des > 0.0)		/* one real root, one complex conjugate pair */
		{
			double d2 = sqrt(des);
			double t1 = -b2/2 + d2;
			double t2 = -b2/2 - d2;
			double a3;
			double b3;
			Zero_test2 (t1, b2/2, d2);
			Zero_test2 (t2, b2/2, d2);
			a3 = cube_root (t1);
			b3 = cube_root (t2);
			z1 = a3 + b3;
			Zero_test2 (z1, a3, b3);
			z2 = - z1/2;
			t1 = a3-b3;
			Zero_test2 (t1, a3, b3);
			z3 = sqrt(3.0) * t1/2;
			*x1 = z1 - p/3;
			Zero_test2 (*x1, z1, p/3);
			*x2 = z2 - p/3;
			Zero_test2 (*x2, z2, p/3);
			*x3 = z3;
			return (-3);
		}
		else if (des < 0.0)	/* three unequal real roots */
		{
			double temp_r, theta, cos_term, sin_term, t1;

			t1 = b2*b2/4 - des;
			Zero_test2 (t1, b2*b2/4, des);
			temp_r = cube_root (sqrt (b2*b2/4 - des));
			theta = atan2 (sqrt(-des), (-b2/2));
			cos_term = temp_r * cos (theta/3);
			sin_term = temp_r * sin (theta/3) * sqrt(3.0);
			z1 = 2 * cos_term;
			z2 = - cos_term - sin_term;
			Zero_test2 (z2, cos_term, sin_term);
			z3 = - cos_term + sin_term;
			Zero_test2 (z3, cos_term, sin_term);

			*x1 = z1 - p/3;
			Zero_test2 (*x1, z1, p/3);
			*x2 = z2 - p/3;
			Zero_test2 (*x2, z2, p/3);
			*x3 = z3 - p/3;
			Zero_test2 (*x3, z3, p/3);
			return (3);
		}
		else			/* cannot happen */
		{
			fprintf (stderr, "impossible descriminant in solve_cubic\n");
			return (0);
		}

	}
}


int solve_quadratic ( float a, float b, float c, float* x1, float* x2 )
{
	/* 
	*         2
	* Solve ax  + bx + c = 0
	*
	* Return values: 0 - no solution 
	*               -1 - infinite number of solutions
	*                1 - one real root, possibly of multiplicity 2.
	*                    returned in *x1
	*		    2 - two distinct real roots of mutiplicity 1
	*                    returned in *x1 and *x2;
	*		   -2 - two complex roots (conjugates),
	*		        real part returned in *x1,
	*			imaginary part returned in *x2
	*/

	double d;

	if (a == 0)
		return(solve_linear (b, c, x1));

	d = b*b - 4*a*c;
	Zero_test2 (d, b*b, 4*a*c);

	if (d == 0) 
	{
		*x1 = *x2 = -b/(2*a);
		return (1);
	}

	if (d > 0)
	{
		d = sqrt(d);
		*x1 = (-b + d)/(2*a);
		Zero_test2 (*x1, -b/2*a, d/2*a);
		*x2 = (-b - d)/(2*a);
		Zero_test2 (*x2, -b/2*a, -d/2*a);
		return (2);
	}

	if (d < 0)
	{
		d = sqrt(-d);
		*x1 = -b/(2*a);
		*x2 =  d/(2*a);
		return (-2);
	}
	return 0;
}

int solve_linear ( float a, float b, float* x )
{
	/* 
	* Solve ax + b = 0
	*
	* Return values: 0 - no solution 
	*                1 - one real solution, returned in *x
	*               -1 - infinite number of solutions
	*/

	if (a == 0) 
	{
		if (b == 0)
		{
			*x = 0;
			return (-1);
		}
		else
			return (0);
	}
	else
	{
		*x = -b/a;
		return (1);
	}
}
