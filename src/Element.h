
//*************************************************************************
//*			OSU Flow Vector Field							
//*
//*	Created:	Han-Wei Shen, Liya Li		
//*			The Ohio State University
//*	Date:		06/2005
//*	Element: vertex, vertexTopo, face, tetra							 
//
//*************************************************************************/


#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "header.h"
#include "VectorMatrix.h"

//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251 4100 4244)

//////////////////////////////////////////////////////////////////////////
// about the advection point
//////////////////////////////////////////////////////////////////////////
typedef struct PointInfo
{
	VECTOR3 phyCoord;
	VECTOR3 interpolant;	// interpolation coefficients
	int	fromCell;	// advecting result from which cell, mainly used for unstructured grid
	int	inCell;		// in which cell

	PointInfo()
	{
		phyCoord.Zero();
		interpolant.Zero();
		fromCell = inCell = -1;
	};

	void Set(VECTOR3& pcoord, VECTOR3& coeff, int fCell, int iCell)
	{
		phyCoord = pcoord;
		interpolant = coeff;
		fromCell = fCell;
		inCell = iCell;
	};

	void Set(PointInfo& pInfo)
	{
		phyCoord = pInfo.phyCoord;
		interpolant = pInfo.interpolant;
		fromCell = pInfo.fromCell;
		inCell = pInfo.inCell;
	};
}PointInfo;

//////////////////////////////////////////////////////////////////////////
// geometry information for vertex
//////////////////////////////////////////////////////////////////////////
class CVertex
{
public:
	VECTOR3 position;		// position

public:
	CVertex(){};
	CVertex(VECTOR3& pos) {position = pos;}
	void SetPos(VECTOR3& pos) {position = pos;}
};

//////////////////////////////////////////////////////////////////////////
// topology information for vertex
//////////////////////////////////////////////////////////////////////////
class TVertex
{
public:
	vector<int> inc_vert;
	vector<int> inc_tetra;
	unsigned char flag;

public:
	TVertex():flag(0){init();}
	~TVertex() {inc_vert.clear(); inc_tetra.clear();}

	inline void init()
	{
		inc_vert.clear();
		inc_tetra.clear();
	}

	inline void remove_vert(int vid)
	{
		for(int i = 0; i < (int)inc_vert.size(); i ++)
		{
			if(vid == inc_vert[i])
			{ 
				inc_vert[i] = inc_vert[inc_vert.size() - 1];
				inc_vert.pop_back();
				break;	
			}
		}
	}

	inline void add_vert(int vid){ 	inc_vert.push_back(vid); }

	inline void add_unique_vert(int vid)
	{
		bool bfind = false;
		for(int i = 0; i < (int)inc_vert.size() && !bfind; i ++)
		{
			if(vid == inc_vert[i]) 
				bfind = true;
		}
		if(!bfind)
			inc_vert.push_back(vid);
	}

	inline int index(int vid)
	{
		for(int i = 0; i < (int)inc_vert.size(); i ++)
		{
			if(vid == inc_vert[i])
				return i;
		}
		return -1;
	}
};

//////////////////////////////////////////////////////////////////////////
// some pre-computable information about tetra
//////////////////////////////////////////////////////////////////////////
typedef struct TetraInfo
{
	float volume;
	MATRIX3 p2nMatrix;	// matrix used to convert physical coordinate to natural coordinate
}TetraInfo;

//////////////////////////////////////////////////////////////////////////
// basic information for face
//////////////////////////////////////////////////////////////////////////
class CTetra
{
public:
	int ver[4];			// four vertices of this tetra
	int tetra[4];			// four neighboring tetras, for example tetra[0]
					// keeps the neighboring tetra of the vertex face 0 

public:
	CTetra()
	{
		for(int i = 0; i < 4; i ++)
		{
			ver[i] = -1;
			tetra[i] = -1;
		}
	}

	CTetra(int v1, int v2, int v3, int v4)
	{
		ver[0] = v1; ver[1] = v2; ver[2] = v3; ver[3] = v3;
		tetra[0] = tetra[1] = tetra[2] = tetra[3] = -1;
	}

	int index(int v)
	{ 
		if(ver[0] == v) return 0;
		else if(ver[1] == v) return 1;
		else if(ver[2] == v) return 2;
		else if(ver[3] == v) return 3;
		else return -1;
	}

	int tindex(int t)
	{ 
		if(tetra[0] == t) return 0;
		else if(tetra[1] == t) return 1;
		else if(tetra[2] == t) return 2;
		else if(tetra[3] == t) return 3;
		else return -1;
	}
};

//////////////////////////////////////////////////////////////////////////
// public functions
////////////////////////////////////////////////////////////////////////
void ConstructTetraVolume(CTetra* pTetra,TVertex* pVertexTopo,
						  const int verNum, const int tetraNum,
						  const int* pData,
						  const bool bVerTopoOn,const bool bTetraTopoOn);

#endif
