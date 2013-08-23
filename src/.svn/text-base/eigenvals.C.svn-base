/*
THIS SOFTWARE BELONGS TO THE AMERICAN PEOPLE.  IT WAS DEVELOPED
AT THEIR EXPENSE AT NASA AMES RESEARCH CENTER.

NAME: Eigenvalues

REQUIREMENTS: Calculate the eigenvalues of a 3x3 matrix.

DISCUSSION:

HISTORY
Modified for TOPOLOGY				Jan 06				Liya Li
Originated as part of eplane		????				Creon Levit
Modified for TOPOLOGY		 		late 90/early 91	Creon Levit
*/

#include "Topology.h"

int compute_eigenvalues( float m[3][3], float eigenvalues[3] )
{
  float determinant =  m[0][0]*m[1][1]*m[2][2] - m[0][2]*m[1][1]*m[2][0]
                      + m[0][1]*m[1][2]*m[2][0] - m[0][1]*m[1][0]*m[2][2]
					  + m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1];

  float minors = -(  m[0][1]*m[1][0] - m[0][0]*m[1][1] + m[0][2]*m[2][0]
		   + m[1][2]*m[2][1] - m[0][0]*m[2][2] - m[1][1]*m[2][2]);

  float trace = m[0][0] + m[1][1] + m[2][2];

  return (solve_cubic (1.0, -trace, minors, -determinant, &eigenvalues[0], &eigenvalues[1], &eigenvalues[2]));
}
