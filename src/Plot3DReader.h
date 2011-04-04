//////////////////////////////////////////////////////////////////////////
// File Reader for Plot3D format
//
// Created:	Han-Wei Shen, Liya Li
// Date:	06/04/2005
//////////////////////////////////////////////////////////////////////////

#ifndef _PLOT3D_READER_H_
#define _PLOT3D_READER_H_

#include "Field.h"
#include "CurvilinearGrid.h" //added by lijie
//Silence an annoying and unnecessary compiler warning
#pragma warning(disable : 4251)

//using namespace OSUFlow;

//////////////////////////////////////////////////////////////////////////
// assume the grid is not changed for time-varying dataset
// 1 grid file, 1 or more solution files
//////////////////////////////////////////////////////////////////////////

class CPlot3DReader
{
public:
	char *m_strSlnFN;
	char *m_strNodeFN;
	char *m_strTetraFN;
	int m_nTimeSteps;
	int m_nFromTime;
	int m_nToTime;
	int m_nTimevarying;
	int m_nNodeNum;
	int m_nTetraNum;

public:
    CPlot3DReader(const char *slnFN, const char *nodeFN, const char*tetraFN);			// only one time step
	CPlot3DReader(const char *fn);									// using .txt file
	~CPlot3DReader();

	// solution related
	Solution* CreateSolution(void);
	VECTOR3** ReadSolution(int& nodeNum);
	void GetDataFromSln(VECTOR3*, const float*, const float*, const float*, const float*, const int);

	// grid relatetd
	void GetSize(int& nodeNum, int& tetraNum);
	IrregularGrid* CreateIrregularGrid(bool bVerTopoOn, bool bTetraTopoOn);
	void ReadCurvilinearGrid(void);
	void ReadIrregularGrid(bool bVerTopoOn, bool bTetraTopoOn, int nodeNum, int tetraNum, CVertex* pVertexGeom, 	CTetra* pTetra, TVertex* pVertexTopo);
	//added by lijie start to handle curvilinear grid
	void GetDims(int* dim);
	CurvilinearGrid* CreateCurvilinearGrid();
	void ReadCurvilinearGrid( CVertex* pVertexGeom);
	//added by lijie end
	
};

#endif
