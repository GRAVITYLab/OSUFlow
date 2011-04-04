/////////////////////////////////////////////////////////////////////////////
//
//                 OSU Flow Vis Library
//                 Created: Han-Wei Shen, Liya Li 
//                 The Ohio State University	
//                 Date:		06/2005
//
///////////////////////////////////////////////////////////////////////////////


#include "Field.h"

#pragma warning(disable : 4251 4100 4244 4101)


//////////////////////////////////////////////////////////////////////////
//
// CVectorField class definition
//
//////////////////////////////////////////////////////////////////////////
CVectorField::CVectorField()
{
	Reset();
}

CVectorField::CVectorField(Grid* pGrid, Solution* pSolution, int timesteps, 
			   int min_t)
{
	assert((pGrid != NULL) && (pSolution != NULL));
	m_pGrid = pGrid;
	m_pSolution = pSolution;
	m_nTimeSteps = timesteps;
	m_bIsNormalized = false;
	m_MinT = min_t;  m_MaxT = min_t + timesteps -1; 
}

CVectorField::~CVectorField()
{
	if(m_pGrid != NULL)	delete m_pGrid;
	if(m_pSolution != NULL) delete m_pSolution;
}

void CVectorField::Reset(void)
{
	m_nTimeSteps = 0;
	m_pGrid = NULL;
	m_pSolution = NULL;
	m_bIsNormalized = false;
	m_MinT = m_MaxT = -1; 
}

//////////////////////////////////////////////////////////////////////////
// whether field is time varying
//////////////////////////////////////////////////////////////////////////
bool CVectorField::isTimeVarying(void)
{
	return (m_nTimeSteps > 1);
}

//////////////////////////////////////////////////////////////////////////
// to obtain one or more field values at cell vertices
//
// input
//		cellId:			from which cell to obtain value
//		eCellTopoType:	the cell type
//		t:				which time step
// output
//		vNodeData:	include one or more node values
// return
//		1:			operation successful
//		-1:			operation fail
//////////////////////////////////////////////////////////////////////////
int CVectorField::at_cell(int cellId,
						  CellTopoType eCellTopoType,
						  const float t,
						  vector<VECTOR3>& vNodeData)
{
	VECTOR3 nodeData;
	vector<int> vVerIds;
	int iFor;

	vNodeData.clear();
	switch(eCellTopoType)
	{
	case T0_CELL:								// cellId is the node Id
		if(m_pSolution->GetValue(cellId, t, nodeData) == -1)
			return -1;
		vNodeData.push_back(nodeData);
		break;

	case T1_CELL:								// cellId is the edge Id
	case T2_CELL:								// cellId is the face id
	case T3_CELL:								// cellId is tetra or cube id
		if(m_pGrid->getCellVertices(cellId, eCellTopoType, vVerIds) == -1)
			return -1;

		for(iFor = 0; iFor < (int)vVerIds.size(); iFor++)
		{
			if(m_pSolution->GetValue(vVerIds[iFor], t, nodeData) == -1)
				return -1;
			vNodeData.push_back(nodeData);
		}
		break;

	default:
		return -1;
		break;
	}

	return 1;
}

//////////////////////////////////////////////////////////////////////////
// to obtain node data at the physical position pos in timestep t
//
// input
//		fromCell:	if not -1, which cell this position is generated from
//		pos:		physical position of node
//		t:			time step
// output
//		nodeData:	data value of this node
// return
//		1:			operation successful
//		-1:			operation fail
//////////////////////////////////////////////////////////////////////////
int CVectorField::at_phys(VECTOR3 pos, float t, VECTOR3& vecData)
{
	vector<VECTOR3> vNodeData;
	PointInfo pInfo;

	// find the cell this position belongs to
	pInfo.Set(pos, pInfo.interpolant, -1, -1);
	if(m_pGrid->phys_to_cell(pInfo) == -1)
		return -1;

	// get vertex value at cell vertices
	if(at_cell(pInfo.inCell, T3_CELL, t, vNodeData) == -1)
		return -1;

	// interpolate in the cell
	m_pGrid->interpolate(vecData, vNodeData, pInfo.interpolant);

	return 1;
}

int CVectorField::at_phys(const int fromCell, 
			  VECTOR3& pos, 
			  PointInfo& pInfo,
			  const float t, VECTOR3& nodeData)
{
	vector<VECTOR3> vNodeData;
		
	// find the cell this position belongs to
	pInfo.Set(pos, pInfo.interpolant, fromCell, -1);
	if(m_pGrid->phys_to_cell(pInfo) == -1) {
		return -1;
	}

	// get vertex value at cell vertices
	if(at_cell(pInfo.inCell, T3_CELL, t, vNodeData) == -1) {
		return -1;
	}

	// interpolate in the cell
	m_pGrid->interpolate(nodeData, vNodeData, pInfo.interpolant);

	return 1;
}

//////////////////////////////////////////////////////////////////////////
// to obtain node data at position (i, j, k) in time t
//
// input
//		(i, j, k):	position in computational space
//		t:			time step
// output
//		dataValue:	data value of this position
// return
//		1:			operation successful
//		-1:			operation fail
//////////////////////////////////////////////////////////////////////////
int CVectorField::at_vert(const int i, 
			  const int j, 
			  const int k, 
			  const float t, 
			  VECTOR3& dataValue)
{
	int xdim, ydim, zdim;
	getDimension(xdim, ydim, zdim);
	return m_pSolution->GetValue((k*ydim*xdim+j*xdim+i), t, dataValue);
}

//////////////////////////////////////////////////////////////////////////
// to obtain node data at one slice aligned with X, Y, or Z direction
//
// input
//		slice:		which slice
//		eSliceType:	aligned with which direction
//		t:			time step
// output
//		vSliceData:	data value of this slice
// return
//		1:			operation successful
//		-1:			operation fail
//////////////////////////////////////////////////////////////////////////
int CVectorField::at_slice(int slice,
			   SliceType eSliceType,
			   const float t,
			   vector<VECTOR3>& vSliceData)
{
	int xdim, ydim, zdim;
	int iFor, jFor, kFor;
	VECTOR3 dataValue;
	int tag;

	getDimension(xdim, ydim, zdim);
	
	// x aligned slice
	if(eSliceType == X_ALIGNED)
	{
		assert((slice >= 0) && (slice <= (xdim-1)));
		iFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
			for(jFor = 0; jFor < ydim; jFor++)
			{
				tag = m_pSolution->GetValue((kFor*ydim*xdim+jFor*xdim+iFor), t, dataValue);
				if(tag == -1)
					return tag;
				vSliceData.push_back(dataValue);
			}
	}

	// y aligned slice
	else if(eSliceType == Y_ALIGNED)
	{
		assert((slice >= 0) && (slice <= (ydim-1)));
		jFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
			for(iFor = 0; iFor < xdim; iFor++)
			{
				tag = m_pSolution->GetValue((kFor*ydim*xdim+jFor*xdim+iFor), t, dataValue);
				if(tag == -1)
					return tag;
				vSliceData.push_back(dataValue);
			}
	}

	// z aligned slice
	else
	{
		assert((slice >= 0) && (slice <= (zdim-1)));
		kFor = slice;
		for(jFor = 0; jFor < ydim; jFor++)
			for(iFor = 0; iFor < xdim; iFor++)
			{
				tag = m_pSolution->GetValue((kFor*ydim*xdim+jFor*xdim+iFor), t, dataValue);
				if(tag == -1)
					return tag;
				vSliceData.push_back(dataValue);
			}
	}

	return 1;
}

//////////////////////////////////////////////////////////////////////////
// to obtain node data at the computational position (i, j, k) in time t
//
// input
//		(i, j, k):	position in computational space
//		t:			time step
// output
//		dataValue:	data value of this position
// return
//		1:			operation successful
//		-1:			operation fail
//////////////////////////////////////////////////////////////////////////
int CVectorField::at_comp(const int i, 
			  const int j, 
			  const int k, 
			  const float t, 
			  VECTOR3& dataValue)
{
	return 1;
}

//////////////////////////////////////////////////////////////////////////
// volume of a cell
//////////////////////////////////////////////////////////////////////////
float CVectorField::volume_of_cell(int cellId)
{
	return m_pGrid->cellVolume(cellId);
}

//////////////////////////////////////////////////////////////////////////
// to normalize the vector field
// input
// bLocal: whether to normalize in each timestep or through all timesteps.
//		   if locally, then divide the magnitude; if globally, then divide
//		   by the maximal magnitude through the whole field
//////////////////////////////////////////////////////////////////////////
bool CVectorField::IsNormalized()
{
	return m_bIsNormalized;
}

void CVectorField::NormalizeField(bool bLocal)
{
	if(IsNormalized() != true)
	{
		m_bIsNormalized = true;
		m_pSolution->Normalize(bLocal);
	}
}

void CVectorField::ScaleField(float scale) 
{
  m_pSolution->Scale(scale); 
} 

//////////////////////////////////////////////////////////////////////////
// to get physical coordinate by interpolation
// input
// cellId: the known cell Id
// eCellTopoType: what kind of cell will be interpolated
// coeff: the interpolation coefficients
// output
// pos: the output physical coordinates
//////////////////////////////////////////////////////////////////////////
int CVectorField::lerp_phys_coord(int cellId, 
				  CellTopoType eCellTopoType, 
				  float* coeff, 
				  VECTOR3& pos)
{
	vector<int> vVertices;
	int iFor;
	VECTOR3 interpCoeff;
	vector<VECTOR3> vData;

	if(m_pGrid->getCellVertices(cellId, eCellTopoType, vVertices) == 0)
		return 0;

	switch(eCellTopoType)
	{
		case T0_CELL:								
			break;
		case T1_CELL:								
			break;
		case T2_CELL:						
			break;
		case T3_CELL:	
			vData.clear();
			for(iFor = 0; iFor < 4; iFor++)
			{
				VECTOR3 v;
				m_pGrid->at_vertex(vVertices[iFor], v);
				vData.push_back(v);
			}

			interpCoeff.Set(coeff[0], coeff[1], coeff[2]);
			m_pGrid->interpolate(pos, vData, interpCoeff);
			break;
		default:
			break;
	}
	return 1;
}

void CVectorField::getDimension(int& xdim, int& ydim, int& zdim)
{
	m_pGrid->GetDimension(xdim, ydim, zdim);
}

void CVectorField::GetInflowRegion(vector<VECTOR3>& inflowVerts, const float t)
{
	VECTOR3 inwardFaceN[6];
	int iFor, jFor, kFor, xdim, ydim, zdim;

	inwardFaceN[0].Set(0,0,1);		// front
	inwardFaceN[1].Set(0,0,-1);		// back
	inwardFaceN[2].Set(0,1,0);		// bottom
	inwardFaceN[3].Set(0,-1,0);		// top
	inwardFaceN[4].Set(1,0,0);		// left
	inwardFaceN[5].Set(-1,0,0);		// right

	getDimension(xdim, ydim, zdim);

	// front
	kFor = 0;
	for(jFor = 0; jFor < ydim; jFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[0]);
			if(dotValue > 0.0)			// inflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				inflowVerts.push_back(index);
			}
		}
	}

	// back
	kFor = zdim-1;
	for(jFor = 0; jFor < ydim; jFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[1]);
			if(dotValue > 0.0)			// inflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				inflowVerts.push_back(index);
			}
		}
	}

	// bottom
	jFor = 0;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[2]);
			if(dotValue > 0.0)			// inflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				inflowVerts.push_back(index);
			}
		}
	}

	// top
	jFor = ydim-1;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[3]);
			if(dotValue > 0.0)			// inflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				inflowVerts.push_back(index);
			}
		}
	}

	// left
	iFor = 0;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(jFor = 0; jFor < ydim; jFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[4]);
			if(dotValue > 0.0)			// inflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				inflowVerts.push_back(index);
			}
		}
	}

	// right
	iFor = xdim-1;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(jFor = 0; jFor < ydim; jFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[5]);
			if(dotValue > 0.0)			// inflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				inflowVerts.push_back(index);
			}
		}
	}
}

void CVectorField::GetOutflowRegion(vector<VECTOR3>& outflowVerts, const float t)
{
	VECTOR3 inwardFaceN[6];
	int iFor, jFor, kFor, xdim, ydim, zdim;

	inwardFaceN[0].Set(0,0,1);		// front
	inwardFaceN[1].Set(0,0,-1);		// back
	inwardFaceN[2].Set(0,1,0);		// bottom
	inwardFaceN[3].Set(0,-1,0);		// top
	inwardFaceN[4].Set(1,0,0);		// left
	inwardFaceN[5].Set(-1,0,0);		// right

	getDimension(xdim, ydim, zdim);

	// front
	kFor = 0;
	for(jFor = 0; jFor < ydim; jFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[0]);
			if(dotValue < 0.0)			// outflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				outflowVerts.push_back(index);
			}
		}
	}

	// back
	kFor = zdim-1;
	for(jFor = 0; jFor < ydim; jFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[1]);
			if(dotValue < 0.0)			// outflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				outflowVerts.push_back(index);
			}
		}
	}

	// bottom
	jFor = 0;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[2]);
			if(dotValue < 0.0)			// outflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				outflowVerts.push_back(index);
			}
		}
	}

	// top
	jFor = ydim-1;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[3]);
			if(dotValue < 0.0)			// outflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				outflowVerts.push_back(index);
			}
		}
	}

	// left
	iFor = 0;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(jFor = 0; jFor < ydim; jFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[4]);
			if(dotValue < 0.0)			// outflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				outflowVerts.push_back(index);
			}
		}
	}

	// right
	iFor = xdim-1;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(jFor = 0; jFor < ydim; jFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[5]);
			if(dotValue < 0.0)			// outflow vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				outflowVerts.push_back(index);
			}
		}
	}
}

void CVectorField::GetTangentialflowRegion(vector<VECTOR3>& tanflowVerts, const float t)
{
	VECTOR3 inwardFaceN[6];
	int iFor, jFor, kFor, xdim, ydim, zdim;

	inwardFaceN[0].Set(0,0,1);		// front
	inwardFaceN[1].Set(0,0,-1);		// back
	inwardFaceN[2].Set(0,1,0);		// bottom
	inwardFaceN[3].Set(0,-1,0);		// top
	inwardFaceN[4].Set(1,0,0);		// left
	inwardFaceN[5].Set(-1,0,0);		// right

	getDimension(xdim, ydim, zdim);

	// front
	kFor = 0;
	for(jFor = 0; jFor < ydim; jFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[0]);
			if(dotValue == 0.0)			// tangent vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				tanflowVerts.push_back(index);
			}
		}
	}

	// back
	kFor = zdim-1;
	for(jFor = 0; jFor < ydim; jFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[1]);
			if(dotValue == 0.0)			// tangent vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				tanflowVerts.push_back(index);
			}
		}
	}

	// bottom
	jFor = 0;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[2]);
			if(dotValue == 0.0)			// tangent vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				tanflowVerts.push_back(index);
			}
		}
	}

	// top
	jFor = ydim-1;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(iFor = 0; iFor < xdim; iFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[3]);
			if(dotValue == 0.0)			// tangent vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				tanflowVerts.push_back(index);
			}
		}
	}

	// left
	iFor = 0;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(jFor = 0; jFor < ydim; jFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[4]);
			if(dotValue == 0.0)			// tangent vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				tanflowVerts.push_back(index);
			}
		}
	}

	// right
	iFor = xdim-1;
	for(kFor = 0; kFor < zdim; kFor++)
	{
		for(jFor = 0; jFor < ydim; jFor++)
		{
			VECTOR3 verValue;
			float dotValue;

			at_vert(iFor, jFor, kFor, t, verValue);
			dotValue = dot(verValue, inwardFaceN[5]);
			if(dotValue == 0.0)			// tangent vertex
			{
				VECTOR3 index;
				index.Set(iFor, jFor, kFor);
				tanflowVerts.push_back(index);
			}
		}
	}
}

void CVectorField::GetInflowSlice(vector<VECTOR3>& inflowVerts,
								  const float t,
								  const int slice,
								  const SliceType eSliceType)
{
	VECTOR3 inwardFaceN[6];
	int iFor, jFor, kFor, xdim, ydim, zdim, index;
	vector<VECTOR3> dataAtSlice;

	inwardFaceN[0].Set(0,0,1);		// front
	inwardFaceN[1].Set(0,0,-1);		// back
	inwardFaceN[2].Set(0,1,0);		// bottom
	inwardFaceN[3].Set(0,-1,0);		// top
	inwardFaceN[4].Set(1,0,0);		// left
	inwardFaceN[5].Set(-1,0,0);		// right

	getDimension(xdim, ydim, zdim);
	at_slice(slice, eSliceType, t, dataAtSlice);

	// x aligned, use left normal
	if(eSliceType == X_ALIGNED)
	{
		iFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
		{
			for(jFor = 0; jFor < ydim; jFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = kFor*ydim + jFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[4]);
				if(dotValue > 0.0)			// inflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					inflowVerts.push_back(index);
				}
			}
		}
	}

	// y aligned, use bottom normal
	else if(eSliceType == Y_ALIGNED)
	{
		jFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
		{
			for(iFor = 0; iFor < xdim; iFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = kFor*xdim + iFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[2]);
				if(dotValue > 0.0)			// inflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					inflowVerts.push_back(index);
				}
			}
		}
	}

	// z aligned, use front normal
	else
	{
		kFor = slice;
		for(jFor = 0; jFor < ydim; jFor++)
		{
			for(iFor = 0; iFor < xdim; iFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = jFor*xdim + iFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[0]);
				if(dotValue > 0.0)			// inflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					inflowVerts.push_back(index);
				}
			}
		}
	}
}

void CVectorField::GetOutflowSlice(vector<VECTOR3>& outflowVerts,
								   const float t,
								   const int slice,
								   const SliceType eSliceType)
{
	VECTOR3 inwardFaceN[6];
	int iFor, jFor, kFor, xdim, ydim, zdim, index;
	vector<VECTOR3> dataAtSlice;

	inwardFaceN[0].Set(0,0,1);		// front
	inwardFaceN[1].Set(0,0,-1);		// back
	inwardFaceN[2].Set(0,1,0);		// bottom
	inwardFaceN[3].Set(0,-1,0);		// top
	inwardFaceN[4].Set(1,0,0);		// left
	inwardFaceN[5].Set(-1,0,0);		// right

	getDimension(xdim, ydim, zdim);
	at_slice(slice, eSliceType, t, dataAtSlice);

	// x aligned, use left normal
	if(eSliceType == X_ALIGNED)
	{
		iFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
		{
			for(jFor = 0; jFor < ydim; jFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = kFor*ydim + jFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[4]);
				if(dotValue < 0.0)			// outflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					outflowVerts.push_back(index);
				}
			}
		}
	}

	// y aligned, use bottom normal
	else if(eSliceType == Y_ALIGNED)
	{
		jFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
		{
			for(iFor = 0; iFor < xdim; iFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = kFor*xdim + iFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[2]);
				if(dotValue < 0.0)			// outflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					outflowVerts.push_back(index);
				}
			}
		}
	}

	// z aligned, use front normal
	else
	{
		kFor = slice;
		for(jFor = 0; jFor < ydim; jFor++)
		{
			for(iFor = 0; iFor < xdim; iFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = jFor*xdim + iFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[0]);
				if(dotValue < 0.0)			// outflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					outflowVerts.push_back(index);
				}
			}
		}
	}
}

void CVectorField::GetTangentialflowSlice(vector<VECTOR3>& tanflowVerts,
										  const float t, 
										  const int slice, 
										  const SliceType eSliceType)
{
	VECTOR3 inwardFaceN[6];
	int iFor, jFor, kFor, xdim, ydim, zdim, index;
	vector<VECTOR3> dataAtSlice;

	inwardFaceN[0].Set(0,0,1);		// front
	inwardFaceN[1].Set(0,0,-1);		// back
	inwardFaceN[2].Set(0,1,0);		// bottom
	inwardFaceN[3].Set(0,-1,0);		// top
	inwardFaceN[4].Set(1,0,0);		// left
	inwardFaceN[5].Set(-1,0,0);		// right

	getDimension(xdim, ydim, zdim);
	at_slice(slice, eSliceType, t, dataAtSlice);

	// x aligned, use left normal
	if(eSliceType == X_ALIGNED)
	{
		iFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
		{
			for(jFor = 0; jFor < ydim; jFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = kFor*ydim + jFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[4]);
				if(dotValue == 0.0)			// tangentflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					tanflowVerts.push_back(index);
				}
			}
		}
	}

	// y aligned, use bottom normal
	else if(eSliceType == Y_ALIGNED)
	{
		jFor = slice;
		for(kFor = 0; kFor < zdim; kFor++)
		{
			for(iFor = 0; iFor < xdim; iFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = kFor*xdim + iFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[2]);
				if(dotValue == 0.0)			// tangentflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					tanflowVerts.push_back(index);
				}
			}
		}
	}

	// z aligned, use front normal
	else
	{
		kFor = slice;
		for(jFor = 0; jFor < ydim; jFor++)
		{
			for(iFor = 0; iFor < xdim; iFor++)
			{
				VECTOR3 verValue;
				float dotValue;

				index = jFor*xdim + iFor;
				verValue = dataAtSlice[index];
				dotValue = dot(verValue, inwardFaceN[0]);
				if(dotValue == 0.0)			// tangentflow vertex
				{
					VECTOR3 index;
					index.Set(iFor, jFor, kFor);
					tanflowVerts.push_back(index);
				}
			}
		}
	}
}

// compute the curl
void CVectorField::curl(float deltaX, 
						float deltaY, 
						float deltaZ,		// grid spacing
						VECTOR3& uL,		// vector value of i-1
						VECTOR3& uR,		// vector value of i+1
						VECTOR3& vL, 
						VECTOR3& vR, 
						VECTOR3& wL, 
						VECTOR3& wR,
						VECTOR3& vort)
{
	float x, y, z;
	float oneoverdeltax, oneoverdeltay, oneoverdeltaz;
	oneoverdeltax = 1.0 / deltaX;
	oneoverdeltay = 1.0 / deltaY;
	oneoverdeltaz = 1.0 / deltaZ;
	
	x = (wR[1] - wL[1])*oneoverdeltay - (vR[2] - vL[2])*oneoverdeltaz;
	y = (uR[2] - uL[2])*oneoverdeltaz - (wR[0] - wL[0])*oneoverdeltax;
	z = (vR[0] - vL[0])*oneoverdeltax - (uR[1] - uL[1])*oneoverdeltay;
	vort.Set(x, y, z);
}

// get the curl at point "pos" at time step t
void CVectorField::at_curl(int t,						// which time step
						   VECTOR3& pos,				// which point
						   VECTOR3& curlAtPos)			// curl generated
{
	float xRatio, yRatio, zRatio;						// interpolation ratio
	VECTOR3 uL, uR, vL, vR, wL, wR;
	bool xBoundary, yBoundary, zBoundary;
	PointInfo pInfo;
	float deltaX, deltaY, deltaZ;
	VECTOR3 p1, p2;
	float gridX, gridY, gridZ;
	int xdim, ydim, zdim;
	
	getDimension(xdim, ydim, zdim);
	m_pGrid->GetGridSpacing(0, gridX, gridY, gridZ);

	xRatio = pos[0] - floor(pos[0]);
	yRatio = pos[1] - floor(pos[1]);
	zRatio = pos[2] - floor(pos[2]);

	if(xRatio == 0.0)					// point on yz plane
	{
		if(pos[0] == 0.0)
		{
			p1.Set(pos[0], pos[1], pos[2]);
			p2.Set(pos[0]+1.0, pos[1], pos[2]);
			deltaX = gridX;
		}
		else
		{
			p1.Set(pos[0]-1.0, pos[1], pos[2]);
			p2.Set(pos[0]+1.0, pos[1], pos[2]);
			deltaX = 2*gridX;
		}
	}
	else if(xRatio == gridX)
	{
		if(pos[0] == (xdim-1))
		{
			p1.Set(pos[0]-1.0, pos[1], pos[2]);
			p2.Set(pos[0], pos[1], pos[2]);
			deltaX = gridX;
		}
		else
		{
			p1.Set(pos[0]-1.0, pos[1], pos[2]);
			p2.Set(pos[0]+1.0, pos[1], pos[2]);
			deltaX = 2*gridX;
		}
	}
	else
	{
		p1.Set(floor(pos[0]), pos[1], pos[2]);
		p2.Set( ceil(pos[0]), pos[1], pos[2]);
		deltaX = gridX;
	}
	at_phys(-1, p1, pInfo, t, uL);
	at_phys(-1, p2, pInfo, t, uR);

	if(yRatio == 0.0)					// point on xz plane
	{
		if(pos[1] == 0.0)
		{
			p1.Set(pos[0], pos[1], pos[2]);
			p2.Set(pos[0], pos[1]+1.0, pos[2]);
			deltaY = gridY;
		}
		else
		{
			p1.Set(pos[0], pos[1]-1.0, pos[2]);
			p2.Set(pos[0], pos[1]+1.0, pos[2]);
			deltaY = 2*gridY;
		}
	}
	else if(yRatio == gridY)
	{
		if(pos[1] == (ydim-1))
		{
			p1.Set(pos[0], pos[1]-1.0, pos[2]);
			p2.Set(pos[0], pos[1], pos[2]);
			deltaY = gridY;
		}
		else
		{
			p1.Set(pos[0], pos[1]-1.0, pos[2]);
			p2.Set(pos[0], pos[1]+1.0, pos[2]);
			deltaY = 2*gridY;
		}
	}
	else
	{
		p1.Set(pos[0], floor(pos[1]), pos[2]);
		p2.Set(pos[0],  ceil(pos[1]), pos[2]);
		deltaY = gridY;
	}
	at_phys(-1, p1, pInfo, t, vL);
	at_phys(-1, p2, pInfo, t, vR);

	if(zRatio == 0.0)					// point on xy plane
	{
		if(pos[2] == 0.0)
		{
			p1.Set(pos[0], pos[1], pos[2]);
			p2.Set(pos[0], pos[1], pos[2]+1.0);
			deltaZ = gridZ;
		}
		else
		{
			p1.Set(pos[0], pos[1], pos[2]-1.0);
			p2.Set(pos[0], pos[1], pos[2]+1.0);
			deltaZ = 2*gridZ;
		}
	}
	else if(zRatio == gridZ)
	{
		if(pos[2] == (zdim-1))
		{
			p1.Set(pos[0], pos[1], pos[2]-1.0);
			p2.Set(pos[0], pos[1], pos[2]);
			deltaZ = gridZ;
		}
		else
		{
			p1.Set(pos[0], pos[1], pos[2]-1.0);
			p2.Set(pos[0], pos[1], pos[2]+1.0);
			deltaZ = 2*gridZ;
		}
	}
	else
	{
		p1.Set(pos[0], pos[1], floor(pos[2]));
		p2.Set(pos[0], pos[1],  ceil(pos[2]));
		deltaZ = gridZ;
	}
	at_phys(-1, p1, pInfo, t, wL);
	at_phys(-1, p2, pInfo, t, wR);

	curl(deltaX, deltaY, deltaZ, uL, uR, vL, vR, wL, wR, curlAtPos);
}

// generate vorticity field
void CVectorField::GenerateVortField(int t,				// which time step
									 bool bToNormalize, // whether to normalize the vort field
									 VECTOR3* pVort		// output buffer
									 )
{
	int iFor, jFor, kFor, index, xdim, ydim, zdim;
	float gridX, gridY, gridZ;
	VECTOR3 uL, uR, vL, vR, wL, wR;
	bool xBoundary, yBoundary, zBoundary;

	getDimension(xdim, ydim, zdim);
	m_pGrid->GetGridSpacing(0, gridX, gridY, gridZ);

	for(kFor = 0; kFor < zdim; kFor++)
		for(jFor = 0; jFor < ydim; jFor++)
			for(iFor = 0; iFor < xdim; iFor++)
			{
				xBoundary = yBoundary = zBoundary = false;
                index = kFor*ydim*xdim+jFor*xdim+iFor;

				if(iFor == 0)
				{
					xBoundary = true;
					at_vert(iFor, jFor, kFor, t, uL);
				}
				else
					at_vert(iFor-1, jFor, kFor, t, uL);

				if(iFor == (xdim-1))
				{
					xBoundary = true;
					at_vert(iFor, jFor, kFor, t, uR);
				}
				else
					at_vert(iFor+1, jFor, kFor, t, uR);

				if(jFor == 0)
				{
					yBoundary = true;
					at_vert(iFor, jFor, kFor, t, vL);
				}
				else
					at_vert(iFor, jFor-1, kFor, t, vL);

				if(jFor == (ydim-1))
				{
					yBoundary = true;
					at_vert(iFor, jFor, kFor, t, vR);
				}
				else
					at_vert(iFor, jFor+1, kFor, t, vR);

				if(kFor == 0)
				{
					zBoundary = true;
					at_vert(iFor, jFor, kFor, t, wL);
				}
				else
					at_vert(iFor, jFor, kFor-1, t, wL);

				if(kFor == (zdim-1))
				{
					zBoundary = true;
					at_vert(iFor, jFor, kFor, t, wR);
				}
				else
					at_vert(iFor, jFor, kFor+1, t, wR);

				curl(xBoundary? gridX: gridX*2.0, 
					 yBoundary? gridY: gridY*2.0,
					 zBoundary? gridZ: gridZ*2.0, 
					 uL, uR, vL, vR, wL, wR, pVort[index]);

				if(bToNormalize)
				{
					float mag = pVort[index].GetMag();
					if(mag > EPS)
					{
						pVort[index][0] /= mag;
						pVort[index][1] /= mag;
						pVort[index][2] /= mag;
					}
				}
			}
}

void CVectorField::GenerateLapField(int t, bool bToNormalize, VECTOR3* pLap)
{
	int iFor, jFor, kFor, index, xdim, ydim, zdim;
	VECTOR3 p1, p2, p3;
	float u, v, w, gridX, gridY, gridZ;
		
	getDimension(xdim, ydim, zdim);
	m_pGrid->GetGridSpacing(0, gridX, gridY, gridZ);
	gridX *= gridX;
	gridY *= gridY;
	gridZ *= gridZ;

	for(kFor = 0; kFor < zdim; kFor++)
		for(jFor = 0; jFor < ydim; jFor++)
			for(iFor = 0; iFor < xdim; iFor++)
			{
				index = kFor*ydim*xdim + jFor*xdim + iFor;
				pLap[index].Set(0.0, 0.0, 0.0);
			}

	for(kFor = 1; kFor < (zdim-1); kFor++)
		for(jFor = 1; jFor < (ydim-1); jFor++)
			for(iFor = 1; iFor < (xdim-1); iFor++)
			{
				index = kFor*ydim*xdim + jFor*xdim + iFor;

				at_vert(iFor+1, jFor, kFor, t, p1);
				at_vert(iFor, jFor, kFor, t, p2);
				at_vert(iFor-1, jFor, kFor, t, p3);
				u = (p1[0] - 2*p2[0] + p3[0])/gridX;
				v = (p1[1] - 2*p2[1] + p3[1])/gridX;
				w = (p1[2] - 2*p2[2] + p3[2])/gridX;
				pLap[index].Set(u, v, w);

				at_vert(iFor, jFor+1, kFor, t, p1);
				at_vert(iFor, jFor, kFor, t, p2);
				at_vert(iFor, jFor-1, kFor, t, p3);
				u = (p1[0] - 2*p2[0] + p3[0])/gridY;
				v = (p1[1] - 2*p2[1] + p3[1])/gridY;
				w = (p1[2] - 2*p2[2] + p3[2])/gridY;
				pLap[index].add(u, v, w);

				at_vert(iFor, jFor, kFor+1, t, p1);
				at_vert(iFor, jFor, kFor, t, p2);
				at_vert(iFor, jFor, kFor-1, t, p3);
				u = (p1[0] - 2*p2[0] + p3[0])/gridZ;
				v = (p1[1] - 2*p2[1] + p3[1])/gridZ;
				w = (p1[2] - 2*p2[2] + p3[2])/gridZ;
				pLap[index].add(u, v, w);
			}
}

