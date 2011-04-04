/*

Created by 
Abon Chaudhuri and Teng-Yok Lee (The Ohio State University)
May, 2010

*/

#include <vector>
using namespace std;

#include "TubeRendererInOpenGL.h"

void 
CTubeRendererInOpenGL::_Draw()
{
	if( 0 == uLid )
		return;

	glPushAttrib(
		GL_POLYGON_BIT |
		0
	);

	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);

//	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, &cMaterial.v4Ambient[0]);
//	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &cMaterial.v4Diffuse[0]);
	glEnable(GL_COLOR_MATERIAL);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, &cMaterial.v4Specular[0]);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, cMaterial.fShininess);

	glColor4fv(&cLine.v4Color[0]);

	switch(iDrawPolygon)
	{
	case DRAW_POLYGON_FILL:
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;

	case DRAW_POLYGON_WIREDFRAME:
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		break;

	default:
		fprintf(stderr, "%s@%s(%d): Unregonized option", __FUNCTION__, __FILE__, __LINE__);
	}
	// MOD-BY-LEETEN 10/01/2010-FROM:
		// glCallList(uLid);
	// TO:
	glPushMatrix(); 
	glCallList(uLid);

	glVertexPointer(4, GL_FLOAT, 0, cVertexArray.pfCoords);
	glNormalPointer(GL_FLOAT, 0, cVertexArray.pfNormals);
	glColorPointer(4, GL_FLOAT, 0, cVertexArray.pfColors);

	glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	if(CColorScheme::COLOR_ON_THE_FLY != this->cColorScheme.iScheme )
	// ADD-BY-LEETEN 01/20/2011-END
		glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	// ADD-BY-LEETEN 01/20/2011-BEGIN
	if(CColorScheme::COLOR_ON_THE_FLY == this->cColorScheme.iScheme )
	{
		for(int t = 0; t < this->vviTracePatchBases.size(); t++)
		{
			bool bIsDrawing;
			_CheckTrace(t, bIsDrawing);
			if( !bIsDrawing )
				continue;
			float fR, fG, fB, fA;
			_GetTraceColor(t, fR, fG, fB, fA);
			glColor4f(fR, fG, fB, fA);
			glDrawElements(GL_TRIANGLES, vviTracePatchLengths[t], GL_UNSIGNED_INT, &cVertexArray.puIndices[vviTracePatchBases[t]]);
		}
	}
	else
	// ADD-BY-LEETEN 01/20/2011-END
		glDrawElements(GL_TRIANGLES, iTotalNrOfPatches, GL_UNSIGNED_INT, cVertexArray.puIndices);
	glPopClientAttrib();	// glPushClientAttrib(GL_CLIENT_VERTEX_ARRAY_BIT);

	glPopMatrix(); 
	// MOD-BY-LEETEN 10/01/2010-END
	glPopAttrib();	// glPushAttrib(GL_POLYGON_BIT);
}

void
CTubeRendererInOpenGL::_Update()
{
	if( !this->pDataSource )
		return;

	CTubeRenderer::_Update();

	if( 0 ==uLid )
		uLid = glGenLists(1);

	glNewList(uLid, GL_COMPILE);

	// ADD-BY-LEETEN 04/15/2010-BEGIN
	float fMaxDim = 0;
	for(int i = 1; i < 3; i++)
		fMaxDim = max(fMaxDim, 
			(cBoundingBox.pv3Corners[1][i] - cBoundingBox.pv3Corners[0][i]) );

	glScalef(
		(cBoundingBox.pv3Corners[1][0] - cBoundingBox.pv3Corners[0][0])/fMaxDim, 
		(cBoundingBox.pv3Corners[1][1] - cBoundingBox.pv3Corners[0][1])/fMaxDim, 
		(cBoundingBox.pv3Corners[1][2] - cBoundingBox.pv3Corners[0][2])/fMaxDim);
	// ADD-BY-LEETEN 04/15/2010-END

	glTranslatef(-1.0f, -1.0f, -1.0f); 
	glScalef(
		2.0f/(cBoundingBox.pv3Corners[1][0] - cBoundingBox.pv3Corners[0][0]), 
		2.0f/(cBoundingBox.pv3Corners[1][1] - cBoundingBox.pv3Corners[0][1]), 
		2.0f/(cBoundingBox.pv3Corners[1][2] - cBoundingBox.pv3Corners[0][2]));
	glTranslatef(-cBoundingBox.pv3Corners[0][0], -cBoundingBox.pv3Corners[0][1], -cBoundingBox.pv3Corners[0][2]); 

	if(cVertexArray.pfCoords) 
		free(cVertexArray.pfCoords);
	cVertexArray.pfCoords = (float*)calloc(sizeof(cVertexArray.pfCoords[0]) * 4, vv3Coords.size());
	for(int p = 0, v = 0; v < vv3Coords.size(); v++)
		for(int i = 0; i < 4; i++, p++)
			if( i < 3 )
				cVertexArray.pfCoords[p] = vv3Coords[v][i];
			else
				cVertexArray.pfCoords[p] = 1.0f;
	glVertexPointer(4, GL_FLOAT, 0, cVertexArray.pfCoords);

	if(cVertexArray.pfNormals) 
		free(cVertexArray.pfNormals);
	cVertexArray.pfNormals = (float*)calloc(sizeof(cVertexArray.pfNormals[0]) * 3, vv3Normals.size());
	for(int p = 0,	v = 0; v < vv3Normals.size(); v++)
		for(int		i = 0; i < 3; i++, p++)
			cVertexArray.pfNormals[p] = vv3Normals[v][i];
	glNormalPointer(GL_FLOAT, 0, cVertexArray.pfNormals);

	// ADD-BY-LEETEN 07/07/2010-BEGIN
	if(cVertexArray.pfColors) 
		free(cVertexArray.pfColors);
	cVertexArray.pfColors = (float*)calloc(sizeof(cVertexArray.pfColors[0]) * 4, vv4Colors.size());
	for(int p = 0,	v = 0; v < vv4Colors.size(); v++)
		for(int		i = 0; i < 4; i++, p++)
			cVertexArray.pfColors[p] = vv4Colors[v][i];
	glColorPointer(4, GL_FLOAT, 0, cVertexArray.pfColors);
	// ADD-BY-LEETEN 07/07/2010-END

	if(cVertexArray.puIndices) 
		free(cVertexArray.puIndices);
	cVertexArray.puIndices = (unsigned int*)calloc(sizeof(cVertexArray.puIndices[0]), viPatchPointIndices.size());
	for(int p = 0; p < viPatchPointIndices.size(); p++)
		cVertexArray.puIndices[p] = (unsigned int)viPatchPointIndices[p];

	// ADD-BY-LEETEN 10/01/2010-BEGIN
	iTotalNrOfPatches = viPatchPointIndices.size();
	// ADD-BY-LEETEN 10/01/2010-END

	viPatchPointIndices.clear();
	vv3Normals.clear();
	vv3Coords.clear();
	// ADD-BY-LEETEN 07/07/2010-BEGIN
	vv4Colors.clear();
	// ADD-BY-LEETEN 07/07/2010-END

	glEndList();
}

CTubeRendererInOpenGL::CTubeRendererInOpenGL(void)
{
}

CTubeRendererInOpenGL::~CTubeRendererInOpenGL(void)
{
}

/*

$Log: TubeRendererInOpenGL.cpp,v $
Revision 1.10  2011/01/20 17:16:30  leeten

[01/19/2010]
1. [DEL] Remove old deleted code segments.
2. [MOD] Enable the color array when CColorScheme::COLOR_ON_THE_FLY != cColorScheme.iScheme.
3. [MOD] Draw the patches per streamlines.

Revision 1.9  2011/01/20 17:15:24  leeten

[01/19/2010]
1. [DEL] Remove old deleted code segments.
2. [MOD] Compute the #points per streamlines in the new _TraverseTraceEnd().
3. [MOD] Enable the color array when CColorScheme::COLOR_ON_THE_FLY != cColorScheme.iScheme.

Revision 1.8  2010/10/01 20:36:09  leeten

[10/01/2010]
1. Checkin the merged version from r.186.

Revision 1.6  2010/08/15 12:34:31  leeten

[08/14/2010]
1. [ADD] Add the authorship to the beginning of the source codes.

Revision 1.5  2010/08/14 23:27:20  leeten

[08/14/2010]
1. [MOD] Change the cast from C++ style to C style.

Revision 1.4  2010/07/07 17:46:12  leeten
[07/07/2010]
1. [ADD] Add a new field vv4Colors to store the vertex color during the traversal of streamlines.

Revision 1.3  2010/07/05 14:27:19  leeten

[07/05/2010]
1. [ADD] Add structures and code segments to supprt the use of vertex array.

Revision 1.2  2010/04/16 17:27:26  leeten

[04/16/2010]
1. [ADD] Scale the object coordinate to keep the aspect ratio.

Revision 1.1.1.1  2010/04/12 21:31:38  leeten

[04/12/2010]
1. [1ST] First Time Checkin.


*/
