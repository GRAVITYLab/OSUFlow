
#include "CurvilinearGrid.h"

const int CurvilinearGrid::subtetra_vertex[10][4] = {
	// even hexahedron cases
	{3, 6, 5, 7}, {0, 3, 5, 1}, {0, 6, 3, 2}, {0, 5, 6, 4}, {0, 5, 3, 6},
		// odd hexahedron cases
	{1, 4, 2, 0}, {1, 2, 7, 3}, {1, 7, 4, 5}, {2, 4, 7, 6}, {1, 2, 4, 7}
};

const int CurvilinearGrid::subtetra_face[10][4][3] = {
	// even hexahedron cases
	{{3, 6, 5}, {3, 5, 7}, {3, 7, 6}, {5, 6, 7}},
	{{0, 3, 5}, {0, 1, 3}, {0, 5, 1}, {1, 5, 3}},
	{{0, 6, 3}, {0, 3, 2}, {0, 2, 6}, {2, 3, 6}},
	{{0, 5, 6}, {0, 4, 5}, {0, 6, 4}, {4, 6, 5}},
	{{0, 5, 3}, {0, 3, 6}, {0, 6, 5}, {3, 5, 6}},
	
	// odd hexahedron cases
	{{1, 4, 2}, {0, 1, 2}, {0, 2, 4}, {0, 4, 1}},
	{{1, 2, 7}, {1, 3, 2}, {1, 7, 3}, {2, 3, 7}},
	{{1, 7, 4}, {1, 4, 5}, {1, 5, 7}, {4, 7, 5}},
	{{2, 4, 7}, {2, 6, 4}, {2, 7, 6}, {4, 6, 7}},
	{{1, 2, 4}, {1, 4, 7}, {1, 7, 2}, {2, 7, 4}}
};

const int CurvilinearGrid::hexa_face[2][6][4] = {
	// even hexahedron case
	{{0, 2, 6, 4}, {3, 1, 5, 7}, {0, 4, 5, 1},
	{3, 7, 6, 2}, {0, 1, 3, 2}, {5, 4, 6, 7}},
	
	// odd hexahedron cases
	{{2, 6, 4, 0}, {1, 5, 7, 3}, {1, 0, 4, 5},
	{2, 3, 7, 6}, {1, 3, 2, 0}, {4, 6, 7, 5}}
};

const Tetra_step CurvilinearGrid::tetra_step[10][4] = {
	// even hexahedron cases
	{ // stepping from tetrahedron 0
		{FEL_DIR_0,     4, 3},
		{FEL_DIR_POS_I, 8, 1},
		{FEL_DIR_POS_J, 7, 1},
		{FEL_DIR_POS_K, 6, 1}
	},
	{ // stepping from tetrahedron 1
		{FEL_DIR_0,     4, 0},
		{FEL_DIR_NEG_K, 7, 3},
		{FEL_DIR_NEG_J, 6, 3},
		{FEL_DIR_POS_I, 5, 2}
	},
	{ // stepping from tetrahedron 2
		{FEL_DIR_0,     4, 1},
		{FEL_DIR_NEG_K, 8, 3},
		{FEL_DIR_NEG_I, 6, 2},
		{FEL_DIR_POS_J, 5, 3}
	},
	{ // stepping from tetrahedron 3
		{FEL_DIR_0,     4, 2},
		{FEL_DIR_NEG_J, 8, 2},
		{FEL_DIR_NEG_I, 7, 2},
		{FEL_DIR_POS_K, 5, 1}
	},
	{ // stepping from tetrahedron 4
		{FEL_DIR_0,     1, 0},
		{FEL_DIR_0,     2, 0},
		{FEL_DIR_0,     3, 0},
		{FEL_DIR_0,     0, 0}
	},
		
		// odd hexahedron cases
	{ // stepping from tetrahedron 5
		{FEL_DIR_0,     9, 0},
		{FEL_DIR_NEG_K, 3, 3},
		{FEL_DIR_NEG_I, 1, 3},
		{FEL_DIR_NEG_J, 2, 3}
	},
	{ // stepping from tetrahedron 6
		{FEL_DIR_0,     9, 2},
		{FEL_DIR_NEG_K, 0, 3},
		{FEL_DIR_POS_I, 2, 2},
		{FEL_DIR_POS_J, 1, 2}
	},
	{ // stepping from tetrahedron 7
		{FEL_DIR_0,     9, 1},
		{FEL_DIR_NEG_J, 0, 2},
		{FEL_DIR_POS_I, 3, 2},
		{FEL_DIR_POS_K, 1, 1}
	},
	{ // stepping from tetrahedron 8
		{FEL_DIR_0,     9, 3},
		{FEL_DIR_NEG_I, 0, 1},
		{FEL_DIR_POS_J, 3, 1},
		{FEL_DIR_POS_K, 2, 1}
	},
	{ // stepping from tetrahedron 9
		{FEL_DIR_0,     5, 0},
		{FEL_DIR_0,     7, 0},
		{FEL_DIR_0,     6, 0},
		{FEL_DIR_0,     8, 0}
	}
};

inline int min(float X, float Y)
{
	return X<Y?X:Y; 
}
inline static int sign(double d)
{
  return d > 0.0 ? 1 : (d < 0.0 ? -1 : 0);
}

int FEL_orient(const VECTOR3& a, const VECTOR3& b,
	       const VECTOR3& c, const VECTOR3& d)
{
  // computed using floats
  //return sign(FEL_dot(d - a, FEL_cross(b - a, c - a)));

  // compute using doubles for fewer round-off problems
  //return sign(FEL_dot(FEL_vector3d(d - a),
  //		      FEL_cross(FEL_vector3d(b - a), FEL_vector3d(c - a))));

  // computed manually for performance, doubles for fewer round-off problems
  double ba_x = b(0) - a(0);
  double ba_y = b(1) - a(1);
  double ba_z = b(2) - a(2);
  double ca_x = c(0) - a(0);
  double ca_y = c(1) - a(1);
  double ca_z = c(2) - a(2);
  return sign((d(0) - a(0)) * (ba_y * ca_z - ca_y * ba_z) +
	      (d(1) - a(1)) * (ca_x * ba_z - ba_x * ca_z) +
	      (d(2) - a(2)) * (ba_x * ca_y - ca_x * ba_y));
}


//sum of ver is odd
int CurvilinearGrid::simplicial_decomposition_odd(VECTOR3& ijk) 
{
	bool res = ((int)(ijk[0] + ijk[1] + ijk[2])) & 1 ? true : false;
	return res;
}


int CurvilinearGrid::up_cells(VECTOR3 input_ijk, Cell& cell,int decompose)
{
	if(decompose==0)//hexahedral
	{
			  for (int i = 0; i < 8; i++) 
			  {
				  VECTOR3 ijk = input_ijk;
				  if (i & 1) 
				  {
					  ijk[0] -= 1;
					  if (ijk[0] < 0) continue;
				  }
				  if (i & 2) 
				  {
					  ijk[1] -= 1;
					  if (ijk[1] < 0) continue;
				  }
				  if (i & 4) 
				  {
					  ijk[2] -= 1;
					  if (ijk[2] < 0) continue;
				  }
				  if (ijk[0] > xdim() - 2) continue;
				  if (ijk[1] > ydim() - 2) continue;
				  if (ijk[2] > zdim() - 2) continue;


				  cell.type=CELL_HEXAHEDRON;
					cell.ijk=ijk;
					cell.subid=-1;

			  }

	}

	int i;
	if (simplicial_decomposition_odd(input_ijk)) 
	{
		//(((i + j + k) & 1) == 0): subid 0~ 4 else 5~9
		static const int octant_tet[8] = {5, 1, 2, 6, 3, 7, 8, 0};//4 and 9 are in the middle
		for (i = 0; i < 8; i++) //check 8 neighbor?
		{
			VECTOR3 ijk = input_ijk;
			if (i & 1) 
			{
				ijk[0] -= 1;
				if (ijk[0] < 0) continue;
			}
			if (i & 2) 
			{
				ijk[1] -= 1;
				if (ijk[1] < 0) continue;
			}
			if (i & 4) 
			{
				ijk[2] -= 1;
				if (ijk[2] < 0) continue;
			}
			if (ijk[0] > xdim() - 2) continue;
			if (ijk[1] > ydim() - 2) continue;
			if (ijk[2] > zdim() - 2) continue;

			cell.type=CELL_TETRAHEDRON;
			cell.ijk=ijk;
			cell.subid=octant_tet[i];
			return 1;//only return one cell
		}
	}
	else {
		static const int octant_tets[8][4] = 
		{
			{1, 2, 3, 4}, {5, 6, 7, 9}, {5, 6, 8, 9}, {0, 1, 2, 4},
			{5, 7, 8, 9}, {0, 1, 3, 4}, {0, 2, 3, 4}, {6, 7, 8, 9}
		};
		for (i = 0; i < 8; i++) 
		{
			VECTOR3 ijk = input_ijk;
			if (i & 1) 
			{
				ijk[0] -= 1;
				if (ijk[0] < 0) continue;
			}
			if (i & 2) 
			{
				ijk[1] -= 1;
				if (ijk[1] < 0) continue;
			}
			if (i & 4) 
			{
				ijk[2] -= 1;
				if (ijk[2] < 0) continue;
			}
			if (ijk[0] > xdim() - 2) continue;
			if (ijk[1] > ydim() - 2) continue;
			if (ijk[2] > zdim() - 2) continue;
				
			cell.type=CELL_TETRAHEDRON;
			cell.ijk=ijk;
			cell.subid=octant_tets[i][0];
			return 1;//only return one cell

		}
	}
	return 1;
	
}//faster than hexahedral

//locate a cell. decompose the cell into 10 tetrahedral and test if the pt is in one of them
//if not, move on to next cell
//a pt is in a tetrahedral if the pt is inside all the 4 faces of a tetrahedral
int CurvilinearGrid::tetrahedral_walk_locate(VECTOR3 phys_pos, Cell prev_cell, Cell& cell)
{
	
	int res;
	bool new_cell = true;
	//cell = prev_cell;//set the first step
	cell.ijk=prev_cell.ijk; cell.subid=prev_cell.subid; cell.type=prev_cell.type;

	int subid = cell.subid;

	if (cell.type == CELL_TETRAHEDRON) 
	{
		cell.type= CELL_HEXAHEDRON;
		cell.subid = -1;
	}
	else 
	{
		subid = !simplicial_decomposition_odd(cell.ijk) ? 4 : 9;
	}
	int face = 0;
	int faces_tested = 0;
	int total_faces_tested = 0;
	int total_faces_tested_threshold = 4 * 5 * (xdim() + ydim() + zdim());
	bool suppressed_step_off_mesh = false;
	bool outside;
	
	VECTOR3 c[8];
	while (faces_tested < 4) 
	{
		total_faces_tested += 1; 
		if (total_faces_tested > total_faces_tested_threshold) 
			return -1;//POINT_LOCATION_STUCK;
		
		if (new_cell) //step to a different cell
		{
			res = coordinates_at_cell(cell.ijk, c);
			if (res != 1) return res;
			new_cell = false;
		}
		
		//test if the point is inside or outside of the texahedra
		int orientation = FEL_orient(c[subtetra_face[subid][face][0]],
		c[subtetra_face[subid][face][1]],
		c[subtetra_face[subid][face][2]],
		phys_pos);
		
		if (left_handed_cells) 
			orientation = -orientation;
		
		
		// FEL_orient can return an orientation of 0 for one of two reasons.
		// First, phys_pos can be coplanar with the face.  Second, the
		// face (and thus the 3-cell) can be degenerate. 
		if (orientation == 0) 
		{
			res = hexahedral_walk_locate(phys_pos, cell, cell);
			return res;
		}
		outside = orientation < 0;
		if (outside) 
		{
			// attempt to take step
			switch(tetra_step[subid][face].dir)
			{
			case FEL_DIR_0:
				break;
			case FEL_DIR_NEG_I:
				if (cell.ijk(0) == 0) 
				{
					suppressed_step_off_mesh = true;
					goto next_subtetra_face;
				}
				cell.ijk[0]-=1;
				new_cell = true;
				break;
			case FEL_DIR_POS_I:
				if (cell.ijk(0) == xdim() - 2) 
				{
					suppressed_step_off_mesh = true;
					goto next_subtetra_face;
				}
				cell.ijk[0] += 1;
				new_cell = true;
				break;
			case FEL_DIR_NEG_J:
				if (cell.ijk(1) == 0) 
				{
					suppressed_step_off_mesh = true;
					goto next_subtetra_face;
				}
				cell.ijk[1] -= 1;
				new_cell = true;
				break;
			case FEL_DIR_POS_J:
				if (cell.ijk(1) == ydim() - 2) 
				{
					suppressed_step_off_mesh = true;
					goto next_subtetra_face;
				}
				cell.ijk[1] += 1;
				new_cell = true;
				break;
			case FEL_DIR_NEG_K:
				if (cell.ijk(2) == 0) 
				{
					suppressed_step_off_mesh = true;
					goto next_subtetra_face;
				}
				cell.ijk[2] -= 1;
				new_cell = true;
				break;
			case FEL_DIR_POS_K:
				if (cell.ijk(2) == zdim() - 2) 
				{
					suppressed_step_off_mesh = true;
					goto next_subtetra_face;
				}
				cell.ijk[2] += 1;
				new_cell = true;
				break;
			default:
				abort();
			}
			// take step
			int prev_subid = subid;
			int prev_face = face;
			subid = tetra_step[prev_subid][prev_face].tetra;
			face = tetra_step[prev_subid][prev_face].face;
			faces_tested = 0;
			suppressed_step_off_mesh = false;
		}
next_subtetra_face:
		face = (face + 1) & 3;
		faces_tested += 1;
	}
	

//	res = suppressed_step_off_mesh ? FEL_POINT_LOCATE_WALKED_OFF_MESH : 1;
	res = suppressed_step_off_mesh ? -1 : 1;
	
done:
	/*
	if (res == 1) 
	{
		int combined_iblank;
		res = combined_iblank_at_cell(*cell, &combined_iblank);
		if (res == 1) 
		{
			if (return_iblank)
				res = combined_iblank;
			else
				res = (combined_iblank & FEL_PLOT3D_HAS_IBLANK_0) ? 0 : 1;
		}
	}
	
	*/
	return res;
	
} 

// Compute orientation of e with respect to quadrilateral abcd
// by treating abcd as two triangles: abc and acd.
//
// Let there be a viewer who sees abcd counter-clockwise (e.g. from
// a position inside a 3-cell).  Return:
//
//    1 if e on same side of quadrilateral abcd as viewer
//   -1 if e on opposite side of quadrilateral abcd with respect to viewer
//    0 if e is coplanar with abcd
//
// The result is computed by doing two orientation tests, one for e
// with respect to abc, and one for e with respect to acd.
// If these two tests agree, we're done, if they don't, then we
// do more work to figure out what case we're in.
//
int FEL_orient(const VECTOR3& a, const VECTOR3& b,
	       const VECTOR3& c, const VECTOR3& d,
	       const VECTOR3& e)
{
  const int verbosity = 0;
  VECTOR3 ba=(b - a);
  VECTOR3 ca=(c - a);
  VECTOR3 da=(d - a);
  VECTOR3 ea=(e - a);

  int abc_orient = sign(dot(ea, cross(ba, ca)));
  int acd_orient = sign(dot(ea, cross(ca, da)));

  if (abc_orient == acd_orient)
    return abc_orient;

  // now the fun begins ...


  // find shortest edge in each triangle
  double ba_length = ba.GetMag();
  double ca_length = ca.GetMag();
  double da_length = da.GetMag();
  VECTOR3 bc=b-c;
  VECTOR3 cd=c-d;
  double bc_length = bc.GetMag();
  double cd_length = cd.GetMag();
  double ea_length = ea.GetMag();
  double shortest_abc_edge_length = min(ba_length, min(ca_length, bc_length));
  double shortest_acd_edge_length = min(ca_length, min(da_length, cd_length));

  //
  // If point e is much further from a triangle than the shortest edge
  // of the triangle, assume we cannot compute a meaningful orientation
  // result in floating-point.
  //
  // Note that it is still possible that a triangle is not usable
  // for orientation even if the following tests conclude that it is
  // (the three vertices can be colinear yet not close together).
  // We ignore that case.
  // 
  double relative_distance_threshold = 1e6;
  bool abc_usable_for_orientation =
    (shortest_abc_edge_length > 0.0) &&
    (ea_length / shortest_abc_edge_length < relative_distance_threshold);
  bool acd_usable_for_orientation =
    (shortest_acd_edge_length > 0.0) &&
    (ea_length / shortest_acd_edge_length < relative_distance_threshold);

 //    if (!abc_usable_for_orientation && !acd_usable_for_orientation) {
 

  if (!abc_usable_for_orientation && !acd_usable_for_orientation)
    return 0;
  if (!abc_usable_for_orientation)
    return acd_orient;
  if (!acd_usable_for_orientation)
    return abc_orient;

 
  int res;

  //
  // determine whether edge ac of abcd is a ridge or a valley with
  // respect to viewer (who sees abcd counter-clockwise); ridge_valley is:
  //   -1 if ridge, i.e. edge ac protrudes towards viewer
  //    1 if valley, i.e. edge ac protrudes away from viewer
  //    0 abcd is planar
  //
  int ridge_valley = sign(dot(da, cross(ba, ca)));
  
  // If ridge_valley == 0, abc and acd are coplanar, so
  // abc_orient and acd_orient should have agreed above.
  // In this case choose triangle with the normal with the biggest
  // magnitude as the triangle to use for orientation.
  if (ridge_valley == 0) 
  {
	  VECTOR3	abc=cross(ba,ca);
	  VECTOR3	acd=cross(ca,da);

    double abc_normal_magnitude = abc.GetMag();
    double acd_normal_magnitude = acd.GetMag();
    if (abc_normal_magnitude > acd_normal_magnitude)
      res = abc_orient;
    else
      res = acd_orient;
    return res;
  }

  // both triangles assumed non-degenerate, quadrilateral is non-planar,
  // do the case analysis
  switch (((abc_orient + 1)<<4) + ((acd_orient + 1)<<2) + (ridge_valley + 1)) {
  //     abc               acd              r_v
  case ((-1 + 1) << 4) + (( 0 + 1) << 2) + (-1 + 1):  res =  0; break;
  case ((-1 + 1) << 4) + (( 0 + 1) << 2) + ( 1 + 1):  res = -1; break;
  case ((-1 + 1) << 4) + (( 1 + 1) << 2) + (-1 + 1):  res =  1; break;
  case ((-1 + 1) << 4) + (( 1 + 1) << 2) + ( 1 + 1):  res = -1; break;
  case (( 0 + 1) << 4) + ((-1 + 1) << 2) + (-1 + 1):  res =  0; break;
  case (( 0 + 1) << 4) + ((-1 + 1) << 2) + ( 1 + 1):  res = -1; break;
  case (( 0 + 1) << 4) + (( 1 + 1) << 2) + (-1 + 1):  res =  1; break;
  case (( 0 + 1) << 4) + (( 1 + 1) << 2) + ( 1 + 1):  res =  0; break;
  case (( 1 + 1) << 4) + ((-1 + 1) << 2) + (-1 + 1):  res =  1; break;
  case (( 1 + 1) << 4) + ((-1 + 1) << 2) + ( 1 + 1):  res = -1; break;
  case (( 1 + 1) << 4) + (( 0 + 1) << 2) + (-1 + 1):  res =  1; break;
  case (( 1 + 1) << 4) + (( 0 + 1) << 2) + ( 1 + 1):  res =  0; break;
  default:
    res = 0;
    abort();
  }
  return res;
}

//////////////////////////////////////////////////////////////////////////
//
//
//	definition of Regular Cartesian Grid Class
//
//
//////////////////////////////////////////////////////////////////////////
// constructor and deconstructor
CurvilinearGrid::CurvilinearGrid():RegularCartesianGrid()
{
	Reset();
	locate_initialization();
}

CurvilinearGrid::CurvilinearGrid(int* dim, CVertex* pVertexGeom):RegularCartesianGrid(dim[0], dim[1], dim[2])
{	
	m_pVertex=pVertexGeom;
	ComputeBBox();
	locate_initialization();
}

CurvilinearGrid::~CurvilinearGrid()
{
}

void CurvilinearGrid::SetBoundary(VECTOR3& minB, VECTOR3& maxB)
{
	m_vMinBound = minB;
	m_vMaxBound = maxB;

}

void CurvilinearGrid::Reset(void)
{

}

void CurvilinearGrid::Boundary(VECTOR3& minB, VECTOR3& maxB)
{
	minB = m_vMinBound;
	maxB = m_vMaxBound;
}

//////////////////////////////////////////////////////////////////////////
// whether the physical point pos is in bounding box
//////////////////////////////////////////////////////////////////////////
bool CurvilinearGrid::isInBBox(VECTOR3& pos)
{
	if( (pos[0] >= m_vMinBound[0]) && (pos[0] <= m_vMaxBound[0]) &&
		(pos[1] >= m_vMinBound[1]) && (pos[1] <= m_vMaxBound[1]) &&
		(pos[2] >= m_vMinBound[2]) && (pos[2] <= m_vMaxBound[2]))
		return true;
	else
		return false;
}

// compute a default boundary 
void CurvilinearGrid::ComputeBBox(void)
{
	VECTOR3 minB, maxB;

	minB.Set(xcelldim(), ycelldim(), zcelldim());
	maxB.Set(0,0,0);

	int j;
	int x,y,z;
	//x bound
	for(j=0; j<2; j++)
	{
		x=0;
		if(j>0)
			x=xcelldim();

		for(int z=0; z<zdim(); z++)
		{
			for(int y=0;y<ydim();y++)
			{
				VECTOR3 v;
				if(1==coordinates_at_vertex(VECTOR3(x,y,z),&v))
				{
					for(int i=0; i<3; i++)
					{
						if(v[i]<minB[i])
							minB[i]=v[i];
						if(v[i]>maxB[i])
							maxB[i]=v[i];
					}
				}

			}
		}
	}
	//y bound
	for(j=0; j<2; j++)
	{
		y=0;
		if(j>0)
			y=ycelldim();

		for(int z=0; z<zdim(); z++)
		{
			for(int x=0;x<xdim();x++)
			{
				VECTOR3 v;
				if(1==coordinates_at_vertex(VECTOR3(x,y,z),&v))
				{
					for(int i=0; i<3; i++)
					{
						if(v[i]<minB[i])
							minB[i]=v[i];
						if(v[i]>maxB[i])
							maxB[i]=v[i];
					}
				}

			}
		}
	}
	//z bound
	for(j=0; j<2; j++)
	{
		z=0;
		if(j>0)
			z=zcelldim();

		for(int y=0; y<ydim(); y++)
		{
			for(int x=0;x<xdim();x++)
			{
				VECTOR3 v;
				if(1==coordinates_at_vertex(VECTOR3(x,y,z),&v))
				{
					for(int i=0; i<3; i++)
					{
						if(v[i]<minB[i])
							minB[i]=v[i];
						if(v[i]>maxB[i])
							maxB[i]=v[i];
					}
				}

			}
		}
	}

	SetBoundary(minB, maxB);
}

//////////////////////////////////////////////////////////////////////////
// for Cartesian grid, this funcion means whether the physical point is 
// in the boundary
//////////////////////////////////////////////////////////////////////////
bool CurvilinearGrid::at_phys(VECTOR3& pos)
{
	// whether in the bounding box
	if(!isInBBox(pos))
		return false;
	
	return true;
}


//////////////////////////////////////////////////////////////////////////
// get the physical coordinate of the vertex
//
// input:
// verIdx: index of vertex
// output:
// pos: physical coordinate of vertex
//////////////////////////////////////////////////////////////////////////
bool CurvilinearGrid::at_vertex(int verIdx, VECTOR3& pos)
{
	int totalVer = xdim() * ydim() * zdim();
	if((verIdx < 0) || (verIdx >= totalVer))
		return false;
	/*
	zidx = verIdx / (xdim() * ydim());
	yidx = verIdx % (xdim() * ydim());
	yidx = verIdx / xdim();
	xidx = verIdx - zidx * xdim() * ydim() - yidx * xdim();
	*/
	float xpos = m_pVertex[verIdx].position[0];
	float ypos = m_pVertex[verIdx].position[1];
	float zpos = m_pVertex[verIdx].position[2];
	
	// pos.Set((float)xidx, (float)yidx, (float)zidx);
	 pos.Set(xpos,ypos,zpos);
	return true;
}


//////////////////////////////////////////////////////////////////////////
// get volume
// input
// cellId:	which cell
// return the volume of this cell
//////////////////////////////////////////////////////////////////////////
float CurvilinearGrid::cellVolume(int cellId)
{
	float volume;
	VECTOR3 c[8],cell;
	
	get_ijk_of_cell(cellId,cell);
	coordinates_at_cell(cell,c);//return 8 velocity value

	volume=dot(c[7]-c[0], cross(c[1],c[3]-c[5]) + 
			 cross(c[2],c[6]-c[3]) +
			 cross(c[4],c[5]-c[6]))/6.;
	return volume;
}


int CurvilinearGrid::central_difference_neighbors(VECTOR3& vc, int dir) 
{
	switch(dir) {
	case DIR_I:
		{
			if (vc[0] == 0)
			{
				if (vc[0] == xdim()-1)
					return -2;
				return -1;
			}
			else if (vc[0] == xdim()-1)
				return 1;
			else
				return 0;
		}
	case DIR_J:
		{
			if (vc[1] == 0)
			{
				if (vc[1] == ydim()-1)
					return -2;
				return -1;
			}
			else if (vc[1] == ydim()-1)
				return 1;
			else
				return 0;
		}
	case DIR_K:
		{
			if (vc[2] == 0)
			{
				if (vc[2] == zdim()-1)
					return -2;
				return -1;
			}
			else if (vc[2] == zdim()-1)
				return 1;
			else
				return 0;
		}
	default:
		abort();
	}
	
	abort();
	return 0;
}

int CurvilinearGrid::get_ijk_of_cell(int cellId, VECTOR3& ijk)
{
	
	if(cellId<0 || (cellId>=(xcelldim())*(ycelldim())*(zcelldim())))
		return -1;
	ijk[2]=cellId/((xcelldim())*(ycelldim()));
	cellId-=ijk[2]*((xcelldim())*(ycelldim()));
	
		  ijk[1]=cellId/((xcelldim()));
		  cellId-=ijk[1]*(xcelldim());
		  
		  ijk[0]=cellId;
		  return 1;
}

int CurvilinearGrid::get_ijk_of_vertex(int vertexId, VECTOR3& ijk)
{
	
	if(vertexId<0 || (vertexId>=xdim()*ydim()*zdim()))
		return -1;
	ijk[2]=vertexId/(xdim()*ydim());
	vertexId-=ijk[2]*(xdim()*ydim());
	
		  ijk[1]=vertexId/xdim();
		  vertexId-=ijk[1]*xdim();
		  
		  ijk[0]=vertexId;
		  return 1;
}



int	CurvilinearGrid::init_jacobian_vectors( VECTOR3* dp)
{
	//initial some values, need the 8 corners data
	xyz_o   = dp[0];
	xyz_r   = dp[1]-dp[0];
	xyz_s   = dp[2]-dp[0];
	xyz_t   = dp[4]-dp[0];
	xyz_rs  = dp[3]-dp[2]-xyz_r;           // these fail if xyz_* are FEL_vector3d
	xyz_rt  = dp[5]-dp[4]-xyz_r;           
	xyz_st  = dp[6]-dp[4]-xyz_s;
	xyz_rst = dp[1]-dp[3]-dp[5]+dp[7]-xyz_st;
	
	return 1;
}


int	CurvilinearGrid::phys_to_comp_coords(const VECTOR3  p,
										 double& r,
										 double& s,
										 double& t
										 )
										 
{
	
	
	VECTOR3 err[6];
	
	// project r,s,t to physical space
	
	VECTOR3 xp = xyz_o 
		+ r*xyz_r + s*xyz_s + t*xyz_t 	
		+ r*s*xyz_rs + r*t*xyz_rt + s*t*xyz_st 
		+ r*s*t*xyz_rst;
	
	
	
	int i;
	for (i=0;i<5;i++) // try a max of 5 Newton-Raphson iterations
	{
		
		// calculate physical space difference vector
		
		VECTOR3 dp = p-xp;
		
		if (dp.GetMag()<1e-6)
		{
			break;
		}
		
		
		// now construct Jacobian at current r,s,t and invert it
		MATRIX3 mi;
		inverse_jacobian_t(mi,r,s,t);
		
		// hit physical space difference vector with J^-1
		
		VECTOR3 dc = mi*dp;
		
		//  cout << "difference (computational) " << dc << endl;
		
		r += dc[0];
		s += dc[1];
		t += dc[2];
		
		// project new r,s,t to physical space
		
		xp = xyz_o 
			+ r*xyz_r + s*xyz_s + t*xyz_t 	
			+ r*s*xyz_rs + r*t*xyz_rt + s*t*xyz_st 
			+ r*s*t*xyz_rst;
		
		
	}
	
	if (i<5) return 1;  
	
	return 0;
}


int CurvilinearGrid::
inverse_jacobian_t(MATRIX3& mi, double r, double s, double t)
{
	VECTOR3 dxyzdr = xyz_r + s*xyz_rs + t*xyz_rt 
		+ s*t*xyz_rst;
	
	VECTOR3 dxyzds = xyz_s + r*xyz_rs + t*xyz_st 
		+ r*t*xyz_rst;
	
	VECTOR3 dxyzdt = xyz_t + r*xyz_rt + s*xyz_st 
		+ r*s*xyz_rst;
	
	MATRIX3 m(dxyzdr,dxyzds,dxyzdt);
	MATRIX3 tm=m.transpose();
	tm.inverse(mi);
	
	return   1;
}

int CurvilinearGrid::coordinates_at_cell(VECTOR3 cell, VECTOR3* dp)
{

	std::vector<int> verIds;
	
	int cellId=(int)cell[0]+((int)cell[1])*(xcelldim())+((int)cell[2])*(xcelldim())*(ycelldim());
	if(!getCellVertices(cellId,T4_CELL,verIds))
		return -1;

	if(verIds.size()!=8)//sth wrong
		return -1;

	for(int i=0; i<verIds.size(); i++)//get the 8 vertices of cell's id
	{
		VECTOR3 ijk;
		get_ijk_of_vertex(verIds[i], ijk);

		coordinates_at_vertex(ijk, &(dp[i]));
	}
	return 1;
}

int CurvilinearGrid::coordinates_at_vertex(VECTOR3 vertex, VECTOR3* dp)
{
	
	int id=(int)vertex[0]+((int)vertex[1])*xdim()+((int)vertex[2])*xdim()*ydim();
	*dp=m_pVertex[id].position;
	return 1;
}
int CurvilinearGrid::phys_to_cell(PointInfo& pInfo)
{
	//seed points
	//VECTOR3 v;
	//coordinates_at_vertex(VECTOR3(3,2,2),&v);	//return 1 velocity value
	

	//check if it is out of bound
	if(!isInBBox(pInfo.phyCoord))
		return -1;

	int fromcellId=pInfo.fromCell;//.inCell;

	Cell fromcell,incell;
	if(pInfo.fromCell==-1)//first point
	{
		if(-1==locate(pInfo.phyCoord, incell))//pInfo.phyCoord is in cell
			return -1;
		//from cell id = incell
		fromcellId=(int)incell.ijk[0]+((int)incell.ijk[1])*(xcelldim())+((int)incell.ijk[2])*(xcelldim())*(ycelldim());

	}
	else
	{
		VECTOR3 fromijk;
		if(-1==get_ijk_of_cell(pInfo.fromCell, fromijk))
			return -1;

		fromcell.ijk=fromijk; fromcell.subid=-1; fromcell.type=CELL_HEXAHEDRON;

		if(-1==locate(pInfo.phyCoord,fromcell,incell))
		return -1;


	}
	//pInfo.fromCell=cell;

	int cellId=(int)incell.ijk[0]+((int)incell.ijk[1])*(xcelldim())+((int)incell.ijk[2])*(xcelldim())*(ycelldim());
	pInfo.fromCell = fromcellId;
	pInfo.inCell = cellId;

	//convert physical position to computational position
	VECTOR3 dp[8];

	if(1!=coordinates_at_cell(incell.ijk,dp))//return 8 velocity value
		return -1;

	init_jacobian_vectors(dp);
	
	double fCoeff[3];
	fCoeff[0]=fCoeff[1]=fCoeff[2]=0.5;

	if(!phys_to_comp_coords(pInfo.phyCoord, fCoeff[0], fCoeff[1], fCoeff[2]))
		return -1;
	
	float  coeff[3];
	for(int i=0; i<3;i++)
		coeff[i]=fCoeff[i];
	
	pInfo.interpolant.Set(coeff[0],coeff[1],coeff[2]);

	return 1;
}

void CurvilinearGrid::locate_initialization()
{
	//	int res;

	//choose a cell in the approximate middle and test for left-hand-cell
	//up_cell()
	left_handed_cells=false;

	VECTOR3 mid_ver=VECTOR3(xdim() / 2, ydim() / 2, zdim() / 2);

	VECTOR3 c[8];
	Cell mid_cell;
	up_cells(mid_ver,mid_cell,0);

	assert(1 == coordinates_at_cell(mid_cell.ijk, c));
	if (FEL_orient(c[0], c[1], c[3], c[2], c[4]) < 0) 
	{
		left_handed_cells = true;
	}

	// choose initial vp's for starting searches for a close vertex
	// choose center of each of the 6 mesh faces as FEL1 does
	n_initial_locations = 6;
	initial_location = new VECTOR3[n_initial_locations];
	int dim0_050 = (xdim() - 1) / 2;
	int dim0m1 = xdim() - 1;
	int dim1_050 = (ydim() - 1) / 2;
	int dim1m1 = ydim() - 1;
	int dim2_050 = (zdim() - 1) / 2;
	int dim2m1 = zdim() - 1;

	int v = 0;
	initial_location[v++]=VECTOR3(dim0_050, dim1_050, 1);
	initial_location[v++]=VECTOR3(dim0_050, dim1_050, dim2m1);
	initial_location[v++]=VECTOR3(dim0_050, 1, dim2_050);
	initial_location[v++]=VECTOR3(dim0_050, dim1m1, dim2_050);
	initial_location[v++]=VECTOR3(1, dim1_050, dim2_050);
	initial_location[v++]=VECTOR3(dim0m1, dim1_050, dim2_050);

}

//locate a cell close to pp
int CurvilinearGrid::locate_close_vertex_cell(VECTOR3 pp, VECTOR3& v)
{
	// FEL_vertex_cell initial_v;
	// initial_v.set_time(pp.get_time());
	
	
	// Find the v closest to pp from the set of initial locations
	VECTOR3 c,initial_v;
	int closest;
	float min_distance = 3.402823466e+38F;
	
	for (int i = 0; i < n_initial_locations; i++) 
	{
		
		initial_v=initial_location[i];//get index of the 8 vertices
		int res = coordinates_at_vertex(initial_v, &c);//get physical positions of each cell
		
		
		float d2 =	(pp[0]-c[0])*(pp[0]-c[0])+
			(pp[1]-c[1])*(pp[1]-c[1])+
			(pp[2]-c[2])*(pp[2]-c[2]);
		if (d2 < min_distance) {
			min_distance = d2;
			closest = i;//select the closest initial position
		}
		
	}
	
	initial_v=initial_location[closest];
	//start walking from the closest initial
	adaptive_vertex_cell_walk(pp, initial_v, &v);

	return 1;
	
}


//pp: the point to approach		start_v: center of the computational domain, starting from here 	v: result
//a fast walking, first large step size than divided by 2 until step size is 1
//converge to a close v not necessarily the nearest
int	CurvilinearGrid::adaptive_vertex_cell_walk(VECTOR3& pp, VECTOR3& start_v,	VECTOR3* v)
{
	VECTOR3 c;
	VECTOR3 current_v = start_v;
	
	int res = coordinates_at_vertex(start_v, &c);//get physical coords at start_v
	if (res != 1) return res;

	int dim[3];
	dim[0]=xdim(); dim[1]=ydim(); dim[2]=zdim(); 
		
	float current_distance =	(pp[0]-c[0])*(pp[0]-c[0])+
		(pp[1]-c[1])*(pp[1]-c[1])+
		(pp[2]-c[2])*(pp[2]-c[2]);

	float stride_size = 0.1f;
	VECTOR3 stride((int) (xdim() * stride_size),
		(int) (ydim() * stride_size),
		(int) (zdim() * stride_size));
	
	for (;;) {
		bool making_progress = true;
		
		while (making_progress) 
		{
			
			float post_step_distance;
			float best_post_step_distance = current_distance;
			VECTOR3 best_post_step_v;
			
			for (int d = 0; d < 3; d++)
			{
				VECTOR3 post_step_v = current_v;
				post_step_v[d] -= stride[d];//walk in one direction
				if (post_step_v[d] >= 0) 
				{
					
					coordinates_at_vertex(post_step_v, &c);
					
					post_step_distance =	(pp[0]-c[0])*(pp[0]-c[0])+
						(pp[1]-c[1])*(pp[1]-c[1])+
						(pp[2]-c[2])*(pp[2]-c[2]);
					
					if (post_step_distance < best_post_step_distance) 
					{
						best_post_step_distance = post_step_distance;
						best_post_step_v = post_step_v;
					}
				}
				post_step_v = current_v;
				post_step_v[d] += stride[d];//walk in the other direction
				
				if (post_step_v[d] < dim[d]) 
				{
					
					coordinates_at_vertex(post_step_v, &c);
					
					
					post_step_distance =	(pp[0]-c[0])* (pp[0]-c[0])+
						(pp[1]-c[1])* (pp[1]-c[1])+
						(pp[2]-c[2])* (pp[2]-c[2]);
					
					if (post_step_distance < best_post_step_distance) 
					{
						best_post_step_distance = post_step_distance;
						best_post_step_v = post_step_v;
					}
				}
			}
			making_progress = best_post_step_distance < current_distance;//is walking closer?
			if (making_progress) {
				current_v = best_post_step_v;//continue walking
 				current_distance = best_post_step_distance;
			}
		}
		if (stride[0] <= 1 && stride[1] <= 1 && stride[2] <= 1) break;//step too small, stop
		//stride /= 2;
		stride[0] /= 2;//reduce walking step and walk from current position
		stride[1] /= 2;//the walk may oscillate nearby destination before converge
		
		stride[2] /= 2;
	}
	*v = current_v;//return current position
	
	return 1;
	
}



int CurvilinearGrid::hexahedral_walk_locate(VECTOR3 phys_pos, Cell prev_cell, Cell& cell)
{
	  
	int res;
	
	bool new_cell = true;
	cell.ijk = prev_cell.ijk; cell.subid = prev_cell.subid; cell.type = prev_cell.type;
	
	// cell->time = phys_pos.get_time();
	if (cell.type == CELL_TETRAHEDRON) 
	{
		cell.type = CELL_HEXAHEDRON;
		cell.subid = -1;
	}  
	int face = 0;
	int faces_tested = 0;
	int total_faces_tested = 0;
	int total_faces_tested_threshold = 2 * (xdim() + ydim() + zdim());
	
	bool suppressed_step_off_mesh = false;
	bool outside;
	
    int even_odd = !simplicial_decomposition_odd(cell.ijk) ? 0 : 1;
	
	while (faces_tested < 6) 
	{
		total_faces_tested += 1; 
		if (total_faces_tested > total_faces_tested_threshold) 
			return -1;//POINT_LOCATION_STUCK;
		
		
		VECTOR3 c[8];
		if (new_cell) 
		{
			res = coordinates_at_cell(cell.ijk, c);
			if (res != 1) return res;
			new_cell = false;
		}
		
		int orientation = FEL_orient(c[hexa_face[even_odd][face][0]],
			c[hexa_face[even_odd][face][1]],
			c[hexa_face[even_odd][face][2]],
			c[hexa_face[even_odd][face][3]],
			phys_pos);
		
		if (left_handed_cells) 
			orientation = -orientation;
		outside = orientation < 0;
		if (outside) 
		{
			// attempt to take step
			switch(face) {
			case 0:
				if (cell.ijk(0) == 0) {
					suppressed_step_off_mesh = true;
					goto next_hexa_face;
				}
				cell.ijk[0] -= 1;
				new_cell = true;
				break;
			case 1:
				if (cell.ijk(0) == xdim() - 2) {
					suppressed_step_off_mesh = true;
					goto next_hexa_face;
				}
				cell.ijk[0] += 1;
				new_cell = true;
				break;
			case 2:
				if (cell.ijk(1) == 0) {
					suppressed_step_off_mesh = true;
					goto next_hexa_face;
				}
				cell.ijk[1] -= 1;
				new_cell = true;
				break;
			case 3:
				if (cell.ijk(1) == ydim() - 2) {
					suppressed_step_off_mesh = true;
					goto next_hexa_face;
				}
				cell.ijk[1] += 1;
				new_cell = true;
				break;
			case 4:
				if (cell.ijk(2) == 0) {
					suppressed_step_off_mesh = true;
					goto next_hexa_face;
				}
				cell.ijk[2] -= 1;
				new_cell = true;
				break;
			case 5:
				if (cell.ijk(2) == zdim() - 2) {
					suppressed_step_off_mesh = true;
					goto next_hexa_face;
				}
				cell.ijk[2] += 1;
				new_cell = true;
				break;
			default:
				abort();
			}
			
			// take step
			even_odd ^= 1;
			face = face ^ 1;
			faces_tested = 0;
			suppressed_step_off_mesh = false;
		}
next_hexa_face:
		face = (face + 1) % 6;
		faces_tested += 1;
	}
	
	
//	res = suppressed_step_off_mesh ? FEL_POINT_LOCATE_WALKED_OFF_MESH : 1;
	
	res = suppressed_step_off_mesh ? -1 : 1;

	return res;
}

int	CurvilinearGrid::walk_on_x_boundary(VECTOR3 pp, int plane_pos,float& min_d2,VECTOR3& min_v)
{
	//dir:FEL_I;

	//iter.set_time(pp.get_time());
	VECTOR3 iter,offset;
	
	offset=VECTOR3(1,0,0);
	
	int x=plane_pos;
	if(x>0)
	{
		offset=VECTOR3(-1,0,0);
	}
	
	for(int z=0; z<zdim(); z++)
	{
		for(int y=0; y<ydim(); y++)
		{
			
			VECTOR3 v,n;
			VECTOR3 iter=VECTOR3(x,y,z);

			coordinates_at_vertex(iter,&v);
			VECTOR3 p_v = pp - v;
			(void) coordinates_at_vertex(iter + offset, &n);
			VECTOR3 n_v = n - v;
			if (dot(p_v, n_v) < 0.0) continue;
			
			float d2 =(v[0]-pp[0])*(v[0]-pp[0])+
				(v[1]-pp[1])*(v[1]-pp[1])+
				(v[2]-pp[2])*(v[2]-pp[2]);
			
			if (d2 < min_d2) 
			{
				min_d2 = d2;
				min_v = iter;
			}
		}
	}
	return 1;
}

int	CurvilinearGrid::walk_on_y_boundary(VECTOR3 pp, int plane_pos,float& min_d2,VECTOR3& min_v)
{
	//dir:FEL_I;

	//iter.set_time(pp.get_time());
	VECTOR3 iter,offset;
	
	offset=VECTOR3(0,1,0);
	
	int y=plane_pos;
	if(y>0)
	{
		offset=VECTOR3(0,-1,0);
	}
	
	for(int z=0; z<zdim(); z++)
	{
		for(int x=0; x<ydim(); x++)
		{
			
			VECTOR3 v,n;
			VECTOR3 iter=VECTOR3(x,y,z);

			coordinates_at_vertex(iter,&v);
			VECTOR3 p_v = pp - v;
			(void) coordinates_at_vertex(iter + offset, &n);
			VECTOR3 n_v = n - v;
			if (dot(p_v, n_v) < 0.0) continue;
			
			float d2 =(v[0]-pp[0])*(v[0]-pp[0])+
				(v[1]-pp[1])*(v[1]-pp[1])+
				(v[2]-pp[2])*(v[2]-pp[2]);
			
			if (d2 < min_d2) 
			{
				min_d2 = d2;
				min_v = iter;
			}
		}
	}
	return 1;
}

int	CurvilinearGrid::walk_on_z_boundary(VECTOR3 pp, int plane_pos,float& min_d2,VECTOR3& min_v)
{
	//dir:FEL_I;

	//iter.set_time(pp.get_time());
	VECTOR3 iter,offset;
	
	offset=VECTOR3(0,0,1);
	
	int z=plane_pos;
	if(z>0)
	{
		offset=VECTOR3(0,0,-1);
	}
	
	for(int y=0; y<ydim(); y++)
	{
		for(int x=0; x<xdim(); x++)
		{
			
			VECTOR3 v,n;
			VECTOR3 iter=VECTOR3(x,y,z);

			coordinates_at_vertex(iter,&v);
			VECTOR3 p_v = pp - v;
			(void) coordinates_at_vertex(iter + offset, &n);
			VECTOR3 n_v = n - v;
			if (dot(p_v, n_v) < 0.0) continue;
			
			float d2 =(v[0]-pp[0])*(v[0]-pp[0])+
				(v[1]-pp[1])*(v[1]-pp[1])+
				(v[2]-pp[2])*(v[2]-pp[2]);
			
			if (d2 < min_d2) 
			{
				min_d2 = d2;
				min_v = iter;
			}
		}
	}
	return 1;
}
//////////////////////////////////////////////////////////////////////////
// whether the point in the physical position is in the cell
//
// input:
// phyCoord:	physical position
// interpolant:	interpolation coefficients
// output:
// pInfo.interpolant: interpolation coefficient
// return:		returns 1 if in cell
//////////////////////////////////////////////////////////////////////////
bool CurvilinearGrid::isInCell(PointInfo& pInfo, const int cellId)
{
	if(!isInBBox(pInfo.phyCoord))
		return false;
	
	float cx, cy, cz; // computatnoial space x, y, and z 
	cx = (pInfo.phyCoord[0]-m_vMinBound[0])/oneOvermappingFactorX; 
	cy = (pInfo.phyCoord[1]-m_vMinBound[1])/oneOvermappingFactorY; 
	cz = (pInfo.phyCoord[2]-m_vMinBound[2])/oneOvermappingFactorZ; 
	
	int xidx, yidx, zidx;
	xidx = (int)floor(cx); 
	yidx = (int)floor(cy); 
	zidx = (int)floor(cz); 
	
	int inCell = zidx * ycelldim() * xcelldim() + yidx * xcelldim() + xidx;
	if(cellId == inCell)
    {
		pInfo.interpolant.Set(cx - (float)xidx, cy - (float)yidx, cz - (float)zidx);
		return true;
    }
	else
		return true;
}

//////////////////////////////////////////////////////////////////////////
// get vertex list of a cell
// input
//		cellId:		cell Id
//		cellType:	cell type
// output
//		vVertices: the vertex lis of the cell
//////////////////////////////////////////////////////////////////////////
//get vertex ids for each cell
int CurvilinearGrid::getCellVertices(int cellId, 
									 CellTopoType cellType, 
									 vector<int>& vVertices)
{
	int totalCell = xcelldim() * ycelldim() * zcelldim();
	int xidx, yidx, zidx, index;
	
	if((cellId < 0) || (cellId >= totalCell))
		return 0;
	
	vVertices.clear();
	zidx = cellId / (xcelldim() * ycelldim());
	yidx = cellId % (xcelldim() * ycelldim());
	yidx = yidx / xcelldim();
	xidx = cellId - zidx * xcelldim() * ycelldim() - yidx * xcelldim();
	
	for(int kFor = 0; kFor < 2; kFor++)
		for(int jFor = 0; jFor < 2; jFor++)
			for(int iFor = 0; iFor < 2; iFor++)
			{
				index = (zidx+kFor) * ydim() * xdim() + (yidx + jFor) * xdim() + (xidx + iFor);
				vVertices.push_back(index);
			}
			return 1;
}

int CurvilinearGrid::locate(VECTOR3 pp, Cell& prev_cell, Cell& cell) 
{
  int res;

  
  res = tetrahedral_walk_locate(pp, prev_cell, cell);

  if (res!=1) 
    res = locate(pp, cell);
  

  return res;
}

int	CurvilinearGrid::locate(VECTOR3 pp, Cell& cell) 
{
	Cell    initial_cell;
	VECTOR3 bb_lo, bb_hi;
	int n_cells;
	int i, res = 0;
	//1. first test the point against boundingbox
	VECTOR3 minBound, maxBound,vc;
	Boundary(minBound, maxBound);	
	

	//1. first try
	if(pp[0]<minBound[0]||pp[1]<minBound[1]||pp[2]<minBound[2]||
		pp[0]>maxBound[0]||pp[1]>maxBound[1]||pp[2]>maxBound[2])
		
		return -1;//out of bound
	
	//debuging
	
//	vc=VECTOR3(10,10,10);
//VECTOR3	 vex;
//	coordinates_at_cell(vc,&vex);//return 8 velocity value

	//2. second try
	if(1!= locate_close_vertex_cell(pp, vc))
		return -1;
	
	//debug
	//VECTOR3 vex[8];
	//coordinates_at_cell(vc,vex);//return 8 velocity value

	up_cells(vc, initial_cell,0);

	res = tetrahedral_walk_locate(pp, initial_cell, cell);
	//if ((return_iblank && (1 <= res && res <= 15)) || res == 1)
	// goto done;
	if(res==1)//done
		return 1;

	//3. third try
		const int MAX_ADAPTIVE_WALK_DESTINATIONS = 8;
		int n_adaptive_walk_destinations = 0;
		VECTOR3 adaptive_walk_destinations[MAX_ADAPTIVE_WALK_DESTINATIONS];
		adaptive_walk_destinations[n_adaptive_walk_destinations++] = vc;//.get_ijk();
		
		for (i = 0; i < n_initial_locations; i++) 
		{
			vc=initial_location[i];
			adaptive_vertex_cell_walk(pp, vc, &vc);
			bool previously_used_destination = false;
			for (int j = 0; j < n_adaptive_walk_destinations; j++) 
			{
				if (vc == adaptive_walk_destinations[j]) 
				{
					previously_used_destination = true;
					break;
				}
			}
			if (previously_used_destination) 
				continue;
			
			//assert(n_adaptive_walk_destinations < MAX_ADAPTIVE_WALK_DESTINATIONS);
			adaptive_walk_destinations[n_adaptive_walk_destinations++] = vc;//.get_ijk();

			//up_cells(vc, 3, 1, &n_cells, &initial_cell);
			up_cells(vc, initial_cell);

			res = tetrahedral_walk_locate(pp, initial_cell, cell);
//			if ((return_iblank && (1 <= res && res <= 15)) || res == 1)
			if(res==1)
				return res;
		}

	// fourth try

		VECTOR3 v, p_v, n, n_v;
		float d2, min_d2 = 3.402823466e+38F;
		VECTOR3 min_v;//cell
		

		walk_on_x_boundary(pp, 0, min_d2, min_v);
		walk_on_x_boundary(pp, xdim()-1, min_d2, min_v);

		walk_on_y_boundary(pp, 0, min_d2, min_v);
		walk_on_y_boundary(pp, ydim()-1, min_d2, min_v);
		

		walk_on_z_boundary(pp, 0, min_d2, min_v);
		walk_on_z_boundary(pp, zdim()-1, min_d2, min_v);

	
		//up_cells(min_v, 3, 1, &n_cells, &initial_cell);
		up_cells(min_v, initial_cell);

		res = tetrahedral_walk_locate(pp, initial_cell, cell);
  
}

