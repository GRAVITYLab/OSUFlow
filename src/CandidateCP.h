#ifndef _CANDIDATE_H_
#define _CANDIDATE_H_

typedef struct
{
	int location[3];
	float vector[2][2][2][3];

	/* these relate to the sub-cube found by bisection */
	/* usually invalid when self pointed to by tCriticalPoint */
	float cubeLocation[2][2][2][3];
	float cubeSize;
}*tCandidate, SizeOfCandidate;

#endif
