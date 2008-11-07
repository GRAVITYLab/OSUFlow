/*************************************************************************
*						OSU Flow Vector Field							 *
*																		 *
*																		 *
*	Created:	Han-Wei Shen, Liya Li									 *
*				The Ohio State University								 *
*	Date:		06/2005													 *
*																		 *
*	Element: vertex, vertexTopo, face, tetra							 *
*************************************************************************/

#include "Element.h"

#pragma warning(disable : 4251 4100 4244 4101)

//////////////////////////////////////////////////////////////////////////
// construct the tetra and vertex topology information
// input
// pData: the array including index of each tetra
// 
// output
// 
//////////////////////////////////////////////////////////////////////////
void ConstructTetraVolume(CTetra* pTetra,
						  TVertex* pVertexTopo,
						  const int verNum,
						  const int tetraNum,
						  const int* pData, 
						  const bool bVerTopoOn,
						  const bool bTetraTopoOn)
{
	int iFor, jFor;
	CTetra *pOneTetra;
	TVertex* pVTopo;

	if(!bVerTopoOn)
		pVTopo = new TVertex[verNum];
	else
		pVTopo = pVertexTopo;
	
	for(iFor = 0; iFor < tetraNum; iFor++)
	{
		pOneTetra = &(pTetra[iFor]);

		pOneTetra->ver[0] = pData[4*iFor+0];			// vertex index of tetra
		pOneTetra->ver[1] = pData[4*iFor+1];
		pOneTetra->ver[2] = pData[4*iFor+2];
		pOneTetra->ver[3] = pData[4*iFor+3];
		
		// construct the vertex topology
		for(jFor = 0; jFor < 4; jFor++)
			pVTopo[pOneTetra->ver[jFor]].inc_tetra.push_back(iFor);

		// vertex incident information
		pVTopo[pOneTetra->ver[0]].add_unique_vert(pOneTetra->ver[1]);
		pVTopo[pOneTetra->ver[0]].add_unique_vert(pOneTetra->ver[2]);
		pVTopo[pOneTetra->ver[0]].add_unique_vert(pOneTetra->ver[3]);
		pVTopo[pOneTetra->ver[1]].add_unique_vert(pOneTetra->ver[0]);
		pVTopo[pOneTetra->ver[1]].add_unique_vert(pOneTetra->ver[2]);
		pVTopo[pOneTetra->ver[1]].add_unique_vert(pOneTetra->ver[3]);
		pVTopo[pOneTetra->ver[2]].add_unique_vert(pOneTetra->ver[0]);
		pVTopo[pOneTetra->ver[2]].add_unique_vert(pOneTetra->ver[1]);
		pVTopo[pOneTetra->ver[2]].add_unique_vert(pOneTetra->ver[3]);
		pVTopo[pOneTetra->ver[3]].add_unique_vert(pOneTetra->ver[0]);
		pVTopo[pOneTetra->ver[3]].add_unique_vert(pOneTetra->ver[1]);
		pVTopo[pOneTetra->ver[3]].add_unique_vert(pOneTetra->ver[2]);
	}

	// iterate all tetras
	for(iFor= 0; iFor < tetraNum; iFor++)
	{
		for(jFor = 0; jFor < 4; jFor++)
		{
			int vid1 = pTetra[iFor].ver[jFor];
			int vid2 = pTetra[iFor].ver[(jFor + 1) % 4];
			int vid3 = pTetra[iFor].ver[(jFor + 2) % 4];

			TVertex& tv = pVTopo[vid1];
			// search the adjacent tetras
			for(vector<int>::iterator fit = tv.inc_tetra.begin(); fit != tv.inc_tetra.end(); fit++)
			{
				if(*fit <= iFor)
					continue;

				int i1, i2, i3;

				if( (i1 = pTetra[*fit].index(vid1)) == -1 )
					continue;
				if( (i2 = pTetra[*fit].index(vid2)) == -1 )
					continue;
				if( (i3 = pTetra[*fit].index(vid3)) == -1 )
					continue;

				if( pTetra[iFor].tetra[(jFor + 3) % 4] >= 0)
					;
				else
					pTetra[iFor].tetra[(jFor + 3) % 4] = *fit;

				for(int k = 0; k < 4; k ++)
				{
					if((k != i1) && (k != i2) && (k != i3))
					{ 
						if(pTetra[*fit].tetra[k] >= 0)
							;
						else
							pTetra[*fit].tetra[k] = iFor;
						break;
					}
				}
			}
		}
	}

	if(!bVerTopoOn)
		delete[] pVTopo;
}
