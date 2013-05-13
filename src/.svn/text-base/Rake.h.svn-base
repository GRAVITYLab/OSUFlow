/////////////////////////////////////////////////////////////////////////////////////////////
//	File:		Rake.h
//
//	Author:		Liya Li
//
//	Date:		July 2005
//
//	Description:	This class is used to generate seed in the flow field. Currently, all rakes
//					are axis-aligned.
/////////////////////////////////////////////////////////////////////////////////////////////

#ifndef	_RAKE_h_
#define	_RAKE_h_

#include "header.h"
#include "VectorMatrix.h"
#include "Interpolator.h"

enum RakeDim
{
    LINE,
	PLANE,
	SOLID
};

class SeedGenerator
{
public:
	SeedGenerator(const float min[3], const float max[3], const size_t numSeeds[3]);
	~SeedGenerator();

	void SetRakeDim(void);
	size_t GetRakeDim(void);
	void GetSeeds(VECTOR3* pSeeds, const bool bRandom);
	void GetSeeds(VECTOR3* pSeeds, vector<VECTOR3>& pPointSets);

private:
	float rakeMin[3], rakeMax[3];		// minimal and maximal positions
	size_t numSeeds[3];				// number of seeds
	int rakeDimension;					// 0, 1, 2, 3
};

class Rake
{
public:
	virtual void GenSeedRandom( const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed) = 0;
	virtual void GenSeedRegular(const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed) = 0;
};

class LineRake : public Rake
{
public:
	LineRake();
	void GenSeedRandom(const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed);
	void GenSeedRegular(const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed);
	~LineRake();
};

class PlaneRake : public Rake
{
public:
	PlaneRake();
	void GenSeedRandom(const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed);
	void GenSeedRegular(const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed);
	~PlaneRake();
};

class SolidRake : public Rake
{
public:
	SolidRake();
	void GenSeedRandom(const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed);
	void GenSeedRegular(const size_t numSeeds[3], const float min[3], const float max[3], VECTOR3* pSeed);
	~SolidRake();
};

#endif
