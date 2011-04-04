/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//                 Streamlines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef _GRID_H_
#define _GRID_H_

#include "header.h"
#include "Element.h"
#include "Interpolator.h"

//Silence an annoying and unnecessary compiler warning
//#pragma warning(disable : 4251 4100 4244)

enum CellType
{
	TRIANGLE,
	CUBE,
	POLYGONE,
	TETRAHEDRON
};

// define the cell type
enum CellTopoType
{
	T0_CELL,					// vertex
	T1_CELL,					// edge
	T2_CELL,					// triangle, quarilateral
	T3_CELL,					// tetrahedra, cube
	T4_CELL						// hetrahedra, added by lijie
};

enum SliceType
{
	X_ALIGNED,
	Y_ALIGNED,
	Z_ALIGNED
};

//////////////////////////////////////////////////////////////////////////
//
// base class for grid
//
//////////////////////////////////////////////////////////////////////////
class Grid
{
public:
        Grid() {}; 
	virtual ~Grid(){}; 
	// get the dimension
	virtual void GetDimension(int& xdim, int& ydim, int& zdim) = 0;
	// physical coordinate of vertex verIdx
	virtual bool at_vertex(int verIdx, VECTOR3& pos) = 0;
	// whether the physical point is in the boundary
	virtual bool at_phys(VECTOR3& pos) = 0;			
	// get vertex list of a cell
	virtual int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices) = 0;
	// get the cell id and also interpolating coefficients for the given physical position
	virtual int phys_to_cell(PointInfo& pInfo) = 0;
	// interpolation
	virtual void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, VECTOR3 coeff) = 0;
	// the volume of cell
	virtual float cellVolume(int cellId) = 0;
	// type of cell
	virtual CellType GetCellType(void) = 0;
	// get min and maximal boundary
	virtual void Boundary(VECTOR3& minB, VECTOR3& maxB) = 0;
	// set bounding box
	virtual void SetBoundary(VECTOR3& minB, VECTOR3& maxB) = 0;
	// get grid spacing in x,y,z dimensions
	virtual void GetGridSpacing(int cellId, float& xspace, float& yspace, float& zspace) = 0;
	// boundary intersection
	virtual void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, 
					  VECTOR3& endP,float* stepSize, float oldStepSize) = 0;

protected:
	// reset parameters
	virtual void Reset(void) = 0;
	// compute bounding box
	virtual void ComputeBBox(void) = 0;
	// whether the point is in the bounding box
	virtual bool isInBBox(VECTOR3& pos) = 0;
	// whether in a cell
	virtual bool isInCell(PointInfo& pInfo, const int cellId) = 0;
};

//////////////////////////////////////////////////////////////////////////
//
// Cartesian Grid (Regular and Irregular)
//
//////////////////////////////////////////////////////////////////////////
class CartesianGrid : public Grid
{
public:
	// constructor and destructor
       CartesianGrid(int xdim, int ydim, int zdim);
	CartesianGrid();
	~CartesianGrid();
	inline void GetDimension(int& xdim, int& ydim, int& zdim)
	       {xdim = m_nDimension[0]; ydim = m_nDimension[1]; zdim = m_nDimension[2];}

	// physical coordinate of vertex verIdx
	virtual bool at_vertex(int verIdx, VECTOR3& pos) =0; 
	// whether the physical point is in the boundary
	virtual bool at_phys(VECTOR3& pos) =0; 
	// get vertex list of a cell
	virtual int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices) =0; 
	// get the cell id and also interpolating coefficients for the given physical position
	virtual int phys_to_cell(PointInfo& pInfo) =0; 
	
	// interpolation
	virtual void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, VECTOR3 coeff) =0; 
	// the volume of cell
	virtual float cellVolume(int cellId) = 0; 
	// type of cell
	virtual CellType GetCellType(void) = 0; 
	// get min and maximal boundary
	virtual void Boundary(VECTOR3& minB, VECTOR3& maxB) = 0; 
	// set bounding box
	virtual void SetBoundary(VECTOR3& minB, VECTOR3& maxB) = 0; 
	// get grid spacing in x,y,z dimensions
	virtual void GetGridSpacing(int cellId, float& xspace, float& yspace, float& zspace) = 0; 
	// boundary intersection
	virtual void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, 
					  VECTOR3& endP,float* stepSize, float oldStepSize) = 0; 


protected:
	// reset parameters
	void Reset(void);
	// dimension related
	inline int xdim(void) { return m_nDimension[0];}
	inline int ydim(void) { return m_nDimension[1];}
	inline int zdim(void) { return m_nDimension[2];}
	inline int xcelldim(void) {return (m_nDimension[0] - 1);}
	inline int ycelldim(void) {return (m_nDimension[1] - 1);}
	inline int zcelldim(void) {return (m_nDimension[2] - 1);}

	int m_nDimension[3];				// dimension
	VECTOR3 m_vMinBound, m_vMaxBound;	// min and maximal boundary
};

//////////////////////////////////////////////////////////////////////////
//
// regular cartesian grid
//
//////////////////////////////////////////////////////////////////////////
// map coordinates in computational space to physical space
#define UCGridPhy2Comp(x, y, f) (((x) - (y))*(f))

class RegularCartesianGrid : public CartesianGrid
{
private:
	float mappingFactorX;				// mapping from physical space to computational space
	float mappingFactorY;
	float mappingFactorZ;
	float oneOvermappingFactorX;
	float oneOvermappingFactorY;
	float oneOvermappingFactorZ;
	float gridSpacing;			        // the minimal grid spacing of all dimensions

public:
	RegularCartesianGrid(int xdim, int ydim, int zdim);
	RegularCartesianGrid();
	~RegularCartesianGrid();
	// physical coordinate of vertex verIdx
	bool at_vertex(int verIdx, VECTOR3& pos);
	// whether the physical point is in the boundary
	bool at_phys(VECTOR3& pos);			
	// get vertex list of a cell
	int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices);
	// get the cell id and also interpolating coefficients for the given physical position
	int phys_to_cell(PointInfo& pInfo);
	// interpolation
	  void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, VECTOR3 coeff); 
	// the volume of cell
	float cellVolume(int cellId);
	// cell type
	CellType GetCellType(void) {return CUBE;}
	// set bounding box
	void SetBoundary(VECTOR3& minB, VECTOR3& maxB);
	// get min and maximal boundary
	void Boundary(VECTOR3& minB, VECTOR3& maxB);
	// get grid spacing in x,y,z dimensions
	void GetGridSpacing(int cellId, float& xspace, float& yspace, float& zspace) 
	{ xspace = oneOvermappingFactorX; yspace = oneOvermappingFactorY; zspace = oneOvermappingFactorZ; }
	void BoundaryIntersection(VECTOR3&, VECTOR3&, VECTOR3&, float*, float);

protected:
	void Reset(void);
	// compute bounding box
	void ComputeBBox(void);
	// whether the point is in the bounding box
	bool isInBBox(VECTOR3& pos);
	// whether in a cell
	bool isInCell(PointInfo& pInfo, const int cellId);
};

/* 
//Comment the following code out since they have not been implemented
// 
//////////////////////////////////////////////////////////////////////////
//
//	irregular cartesian grid
//
//////////////////////////////////////////////////////////////////////////
class IrregularCartesianGrid : public CartesianGrid
{
private:
	float* m_pXSpacing;			// space array for x, y, z dimension
	float* m_pYSpacing;
	float* m_pZSpacing;

public:
	IrregularCartesianGrid(int xdim, int ydim, int zdim);
	IrregularCartesianGrid();
	~IrregularCartesianGrid();

protected:
	void Reset(void);
};

//////////////////////////////////////////////////////////////////////////
//
// curvilinear grid
//
//////////////////////////////////////////////////////////////////////////
class CurvilinearGrid : public Grid
{
private:
	int m_nDimension[3];				// dimension

public:
	// constructor and deconstructor
	CurvilinearGrid(int xdim, int ydim, int zdim);
	CurvilinearGrid();
	~CurvilinearGrid();
};

*/ 

//////////////////////////////////////////////////////////////////////////
//
// irregular grid
//
//////////////////////////////////////////////////////////////////////////
class IrregularGrid : public Grid
{
private:
	int m_nNodeNum;						// number of nodes
	int m_nTetraNum;					// number of tetras
	CVertex* m_pVertexGeom;				// geometry of all vertices
	CTetra* m_pTetra;					// tetra
	TetraInfo* m_pTetraInfo;			// pre-computed tetra information
	TVertex* m_pVertexTopo;				// vertex topology
	VECTOR3 m_vMinBound, m_vMaxBound;	// min and maximal boundary
	bool m_bTetraInfoInit;				// whether the tetra information is pre-computed

public:
	// constructor and deconstructor
	IrregularGrid();
	IrregularGrid(int nodeNum, int tetraNum, CVertex* pVertexGeom, CTetra* pTetra, TVertex* pVertexTopo);
	~IrregularGrid();

	// from virtual functions
	void Reset(void);
	void GetDimension(int& xdim, int& ydim, int& zdim) {xdim = m_nNodeNum; ydim = m_nTetraNum; zdim = 0;}
	bool at_vertex(int verIdx, VECTOR3& pos);
	bool at_phys(VECTOR3& pos);
	int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices);
	int phys_to_cell(PointInfo& pInfo);
	void interpolate(VECTOR3& nodeData, vector<VECTOR3>& vData, VECTOR3 coeff);
	float cellVolume(int cellId);
	bool isInCell(PointInfo& pInfo, const int cellId);
	CellType GetCellType(void) {return TETRAHEDRON;}

	void ComputeBBox(void);
	void SetBoundary(VECTOR3& minB, VECTOR3& maxB);
	void Boundary(VECTOR3& minB, VECTOR3& maxB);
	bool isInBBox(VECTOR3& pos);

	// irregular specific functions
	void SetTetraInfoInit(bool bInit);
	bool GetTetraInfoInit(void);
	int nextTetra(PointInfo& pInfo, int tetraId);
	void PreGetP2NMatrix(MATRIX3& m, int cellId);
	bool Physical2NaturalCoord(VECTOR3& nCoord, VECTOR3& pCoord, int cellId);

	void GetGridSpacing(int cellId, float& xspace, float& yspace, float& zspace) {}; 
	virtual void BoundaryIntersection(VECTOR3& intersectP, VECTOR3& startP, 
					  VECTOR3& endP,float* stepSize, float oldStepSize){}; 

};

float getStepSize(VECTOR3& p, VECTOR3& p1, VECTOR3& p2, float oldStepSize);
#endif
