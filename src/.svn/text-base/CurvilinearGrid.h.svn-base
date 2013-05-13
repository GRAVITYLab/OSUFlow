#ifndef _CURVILINEARGRID_H
#define _CURVILINEARGRID_H

#include "Grid.h"

#define DIR_I 0
#define DIR_J 1
#define DIR_K 2

typedef struct {
	int dir;
	int tetra;
	int face;
} Tetra_step;

#define CELL_TETRAHEDRON 1
#define CELL_HEXAHEDRON  2

typedef struct {
	int type;
	VECTOR3  ijk;
	int      subid;
} Cell;

/*

phys_to_cell
  |
  ---|
  |  locate
  |  |
  |   ---locate_close_vertex_cell
  |  |           |
  |  |           ---adaptive_vertex_cell_walk
  |  |
  |  ---tetrahedral_walk_locate
  |              |
  |               ---hexahedral_walk_locate
  |
  |
  ----phys_to_comp_coords


  */
class CurvilinearGrid : public RegularCartesianGrid
{

private:
	float mappingFactorX;				// mapping from physical space to computational space
	float mappingFactorY;
	float mappingFactorZ;
	float oneOvermappingFactorX;
	float oneOvermappingFactorY;
	float oneOvermappingFactorZ;
	float gridSpacing;			        // the minimal grid spacing of all dimensions
	
	VECTOR3 xyz_o,xyz_r,xyz_s,xyz_t,xyz_rs,xyz_rt,xyz_st,xyz_rst;
	VECTOR3* initial_location;
	CVertex* m_pVertex;
	bool left_handed_cells;

public:

	CurvilinearGrid(int xdim, int ydim, int zdim);
	CurvilinearGrid(int* dim, CVertex* pVertexGeom);

	CurvilinearGrid();
	~CurvilinearGrid();
	void Reset();
	// physical coordinate of vertex verIdx
	bool at_vertex(int verIdx, VECTOR3& pos);
	// whether the physical point is on the boundary
	bool at_phys(VECTOR3& pos);			
	// get vertex list of a cell
	int getCellVertices(int cellId, CellTopoType cellType, vector<int>& vVertices);
	// get the cell id and also interpolating coefficients for the given physical position
	int phys_to_cell(PointInfo& pInfo);
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
	//void BoundaryIntersection(VECTOR3&, VECTOR3&, VECTOR3&, float*, float);
	bool isInBBox(VECTOR3& pos);
	void ComputeBBox();
	bool isInCell(PointInfo& pInfo, const int cellId);

	int get_ijk_of_vertex(int vertexId, VECTOR3& ijk);
	int get_ijk_of_cell(int cellId, VECTOR3& ijk);

	/**********************************/
	//re-implemented functions
	int central_difference_neighbors(VECTOR3& vc, int dir);
	
	//interpolate
	int	phys_to_comp_coords(const VECTOR3  p,double& r,	double& s,double& t);
		
	//	int get_ijk_of_vertex(int vertexId, VECTOR3& ijk);
	int init_jacobian_vectors( VECTOR3* dp);
	//  int	jacobian_at_vertex(VECTOR3, MATRIX3*);//not used
	//	int	jacobian_at_vertex_id(int verIdx, MATRIX3& jacobian);//not used
	
	//int inverse_jacobian(MATRIX3& mi, double r, double s, double t);
	int inverse_jacobian_t(MATRIX3& mi, double r, double s, double t);
	
	//point location
	void locate_initialization();
	int locate_close_vertex_cell(VECTOR3 pp, VECTOR3& vc);
	
	int	adaptive_vertex_cell_walk(VECTOR3& pp, VECTOR3& start_v,VECTOR3* v);

	int tetrahedral_walk_locate(VECTOR3 , Cell , Cell&);

	int hexahedral_walk_locate(VECTOR3 , Cell , Cell&);
	
	int	locate(VECTOR3, Cell&) ;
	int locate(VECTOR3 pp, Cell& prev_cell, Cell& cell);
//	int	locate(VECTOR3&, VECTOR3&, VECTOR3*) const;
	
	//auxiliary
	int coordinates_at_cell(VECTOR3 cell,VECTOR3* v);//return 8 velocity value
	int coordinates_at_vertex(VECTOR3 pos,VECTOR3* v);	//return 1 velocity value
	int	walk_on_x_boundary(VECTOR3 pp, int plane_pos,float& min_d2,VECTOR3& min_v);
	int	walk_on_y_boundary(VECTOR3 pp, int plane_pos,float& min_d2,VECTOR3& min_v);
	int	walk_on_z_boundary(VECTOR3 pp, int plane_pos,float& min_d2,VECTOR3& min_v);
	int simplicial_decomposition_odd(VECTOR3& ijk) ;
	int up_cells(VECTOR3 input_ijk, Cell& cell,int decompose=1);

	
	static const int hexa_face[2][6][4];
	static const int subtetra_vertex[10][4];
	static const int subtetra_face[10][4][3];
	static const int quadrilateral_vertex[3][4];
	static const int subtriangle_vertex[20][3];
	static const int edge_vertex[21][2];
	static const Tetra_step tetra_step[10][4];
	int n_initial_locations;
};

#define FEL_DIR_0	0
#define FEL_DIR_NEG_I	1
#define FEL_DIR_POS_I	2
#define FEL_DIR_NEG_J	3
#define FEL_DIR_POS_J	4
#define FEL_DIR_NEG_K	5
#define FEL_DIR_POS_K	6

#define FEL_DIR_I       7
#define FEL_DIR_J       8
#define FEL_DIR_K       9

// The canonical vertex numbering for a structured hexahedron.
// The vertex indices in the tables below are in terms of this
// hexahedron numbering.
//
//		    6________7
//		   /|       /|  ^
//		  / |      / |  |
//		4/_______5/  |  k
//		|  2|___ |___|3
//		|  /     |  /  ^
//		| /      | /  /
//		|/_______|/  j
//		0        1
//		   i ->


#endif
