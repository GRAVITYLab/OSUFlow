#ifndef _NEWTON_H_
#define _NEWTON_H_

#include "header.h"
#include "CandidateCP.h"
#include "Interpolator.h"

typedef struct
{
	int maxIterations;
	float outerLimit;
	float epsilonFactor;
	float minFloatFactor;
	int numberOfBisections;
	int pureBisection;
	float isInCellFuzz;
}*tNewton, SizeOfNewton;

extern tNewton theNewton;

int newton ( tCandidate, tCandidate, float, float, float );
float CELLtrilinear_interp( float*, float, float, float);

#endif
