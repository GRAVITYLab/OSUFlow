#include "PointRenderer.h"

void
CPointRenderer::_TraversePointsBegin()
{
	glPushMatrix(); 

	glTranslatef(-1.0f, -1.0f, -1.0f); 
	glScalef(
		2.0f/(cBoundingBox.pv3Corners[1][0] - cBoundingBox.pv3Corners[0][0]), 
		2.0f/(cBoundingBox.pv3Corners[1][1] - cBoundingBox.pv3Corners[0][1]), 
		2.0f/(cBoundingBox.pv3Corners[1][2] - cBoundingBox.pv3Corners[0][2]));
	glTranslatef(-cBoundingBox.pv3Corners[0][0], -cBoundingBox.pv3Corners[0][1], -cBoundingBox.pv3Corners[0][2]);

	_InitHeadlight();
}

void
CPointRenderer::_InitHeadlight()
{
	GLfloat light_position0[] = { 0, 0, 10 };
	GLfloat light_direction[] = { 0, 0, -1, 0 };
	GLfloat light_ambient[] = { 0, 0, 0, 0 };
	GLfloat light_diffuse[] = { 0.3, 0.3, 0.3, 0 };
	GLfloat light_specular[] = { 0.2, 0.2, 0.2, 0 };
	GLfloat light_expo[] = { 4 };

	glMatrixMode( GL_MODELVIEW );
	glPushMatrix();
	glLoadIdentity();
		
	glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
	//glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, light_direction);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);

	glPopMatrix();

}// end function


void
CPointRenderer::_TraversePointsEnd()
{
	glPopMatrix();
}// end function

void CPointRenderer::_TraversePoints_Point()
{
	printf( "Rendering Points...\n" );
	
	// Initialize local variables
	sl_list = (const list<vtListSeedTrace*>*)this->pDataSource;
	int iT=0;

	// Get total number of particles
	// Also, maximum lifetime of a particle
	maxT = 0;
	particlecount = _CountParticles( &maxT, &numT );	
	printf("Total number of particles: %d\n", particlecount); 
	printf("Max Lifetime of a particle: %d\n", maxT ); 
	printf("Number of traces: %d\n", numT ); 

	// Allocate memory for rendering quads
	if( pointArray != NULL ) delete [] pointArray;
	pointArray = new float[particlecount * 3];
	colorArray = new float[particlecount * 3];

	// Apply transformations
	_TraversePointsBegin();

	// Scan through all points and asssign bins to them
	iT = 0;
	int offset = 0;
	for(list<vtListSeedTrace*>::const_iterator
			pIter =  sl_list->begin(); 
		pIter!=sl_list->end(); 
		pIter++, iT++) 
	{
	    // Get next trace
		const vtListSeedTrace *trace = *pIter; 

		int iP = 0;
		VECTOR3 p, nextp;
		for(list<VECTOR3*>::const_iterator
				pnIter = trace->begin(); 
			pnIter!= trace->end(); 
			pnIter++, iP++) 
		{
			// Get next particle
			p = **pnIter; 

			// Get World Space Coordinates of the particle
			VECTOR4 p4( p );
			p4 = modelview * p4; 

			// Store vertex location
			pointArray[offset*3] = p(0);	pointArray[offset*3+1] = p(1); pointArray[offset*3+2] = p(2);			
			
			// Assign color based on time
			float f = (float)iP / (float)maxT;
			colorArray[offset*3] = f;	colorArray[offset*3+1] = 0.5f; colorArray[offset*3+2] = 1-f;

			offset ++;

		}// end inner for: for each trace

	}// end outer for

	_TraversePointsEnd();

}

void
CPointRenderer::_RenderPoints_Point()
{
	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glEnable (GL_POINT_SMOOTH);

	glClearColor (0.0, 0.0, 0.0, 0.0);
	glPointSize( 2.0f );

	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_COLOR_ARRAY );

	glVertexPointer( 3, GL_FLOAT, 0, pointArray );
	glColorPointer( 3, GL_FLOAT, 0, colorArray );

	if( singletracemodeEnabled == 1)
	{
		// Render the vertex array	
		glDrawArrays( GL_POINTS, currentTraceOffset, particlecountCurrentTrace ) ;
	}
	else
	{
		// Render the vertex array	
		glDrawArrays( GL_POINTS, 0, particlecount ) ;

	}

	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_COLOR_ARRAY );

}// end function

void
CPointRenderer::_ZMinMaxEyeSpace( float *z )
{
	// Convert each object space vertex to eye space to get eye space z-value 
	VECTOR4 v_obj, v_eye;
	float zmin = 10000, zmax = -10000;
	
	// Point 1
	v_obj.Set( cBoundingBox.pv3Corners[0][0], cBoundingBox.pv3Corners[0][1], cBoundingBox.pv3Corners[0][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	// Point 2
	v_obj.Set( cBoundingBox.pv3Corners[1][0], cBoundingBox.pv3Corners[0][1], cBoundingBox.pv3Corners[0][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	// Point 3
	v_obj.Set( cBoundingBox.pv3Corners[1][0], cBoundingBox.pv3Corners[1][1], cBoundingBox.pv3Corners[0][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	// Point 4
	v_obj.Set( cBoundingBox.pv3Corners[0][1], cBoundingBox.pv3Corners[1][1], cBoundingBox.pv3Corners[0][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	// Point 5
	v_obj.Set( cBoundingBox.pv3Corners[0][0], cBoundingBox.pv3Corners[0][1], cBoundingBox.pv3Corners[1][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	// Point 6
	v_obj.Set( cBoundingBox.pv3Corners[1][0], cBoundingBox.pv3Corners[0][1], cBoundingBox.pv3Corners[1][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	// Point 7
	v_obj.Set( cBoundingBox.pv3Corners[1][0], cBoundingBox.pv3Corners[1][1], cBoundingBox.pv3Corners[1][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	// Point 8
	v_obj.Set( cBoundingBox.pv3Corners[0][0], cBoundingBox.pv3Corners[1][1], cBoundingBox.pv3Corners[1][2], 1 );
	v_eye = modelview * (v_obj);

	if( v_eye(2) < zmin ) zmin = v_eye(2);
	if( v_eye(2) > zmax ) zmax = v_eye(2);

	*z = zmin;
	*(z+1) = zmax;

	//delete v_obj;
	//delete v_eye;

}// end function

int CPointRenderer::_CountParticles( int *maxT, int *numT )
{
	const list<vtListSeedTrace*>* sl_list = (const list<vtListSeedTrace*>*)this->pDataSource;
	int particlecount = 0;
	int iT = 0;
	*maxT = 0;
	*numT = 0;
	
	for(list<vtListSeedTrace*>::const_iterator
		pIter =  sl_list->begin(); 
		pIter!=sl_list->end(); 
		pIter++, iT++) 
	{
	    const vtListSeedTrace *trace = *pIter; 
		particlecount += trace->size();
		if( trace->size() > *maxT ) *maxT = trace->size();

		(*numT) ++;

	}// end for

	return particlecount;
}// end function

void CPointRenderer::_CurrentModelview()
{
	double *mv = new double[16];
	
	glMatrixMode( GL_MODELVIEW );
	glGetDoublev( GL_MODELVIEW_MATRIX, mv );
	
	VECTOR4 *v0 = new VECTOR4( mv[0], mv[1], mv[2], mv[3] );
	VECTOR4 *v1 = new VECTOR4( mv[4], mv[5], mv[6], mv[7] );
	VECTOR4 *v2 = new VECTOR4( mv[8], mv[9], mv[10], mv[11] );
	VECTOR4 *v3 = new VECTOR4( mv[12], mv[13], mv[14], mv[15] );
	
	modelview[0] = *v0;
	modelview[0] = *v1;
	modelview[0] = *v2;
	modelview[0] = *v3;
	
	delete [] mv;
	delete v0, v1, v2, v3;

}// end function

void CPointRenderer::_Draw()
{
	if( 0 == uPid )
		return;

	glCallList(uPid);
}

void CPointRenderer::_Update()
{
	if( NULL == pDataSource )
		return;

	if( 0 ==uPid )
		uPid = glGenLists(1);

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	_TraversePoints_Point();
	_RenderPoints_Point();

	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();
}

void CPointRenderer::_UpdateTraceMode()
{
	if( NULL == pDataSource )
		return;

	if( singletracemodeEnabled == 0 )	singletracemodeEnabled = 1;
	else singletracemodeEnabled = 0;

	if( 0 ==uPid )
		uPid = glGenLists(1);

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	_RenderPoints_Point();

	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();

}

void CPointRenderer::_UpdateSelectedTrace()
{
	if( NULL == pDataSource )
		return;

	if( singletracemodeEnabled == 0 )
		return;

	// Update Current Trace
	currentTrace = ( currentTrace+1 ) % numT;

	list<vtListSeedTrace*>::const_iterator pIter =  sl_list->begin(); 
	currentTraceOffset = 0;

	// Number of particle up to this trace
	for( int i=0; i<currentTrace; i++ )
	{
		currentTraceOffset += (*pIter)->size();
		pIter ++;
	}
		
	// Number of particles in the current trace
	//pIter++;
	particlecountCurrentTrace = (*pIter)->size();

	if( 0 ==uPid )
		uPid = glGenLists(1);

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	_RenderPoints_Point();

	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();

	printf( "%d, %d, %d\n", currentTrace, particlecountCurrentTrace, currentTraceOffset );

}

CPointRenderer::CPointRenderer(void)
{
	// Initialize data structure
	pointArray = NULL;
	colorArray = NULL;
	normalArray = NULL;	

	singletracemodeEnabled = 0;
	currentTrace = -1;
	currentTraceOffset = 0;
	numT = 0;
}

CPointRenderer::~CPointRenderer(void)
{
}
