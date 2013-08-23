#include "PointRendererGlyph.h"

void
CPointRendererGlyph::_TraversePoints_Glyph_Arrow()
{
	printf( "Rendering Arrow-like Glyphs with lighting...\n" );
	
	// Initialize local variables
	const list<vtListSeedTrace*>* sl_list = (const list<vtListSeedTrace*>*)this->pDataSource;
	int iT=0;
	int bincount = 20;
	float *z = new float[2];
	int *offsetToBins = new int[bincount];
	int *pcountOfBins = new int[bincount];
	for(int i=0; i<bincount; i++)
	{
		offsetToBins[i] = 0;
		pcountOfBins[i] = 0;
	}

	// Get total number of particles to begin
	maxT = 0;
	particlecount = _CountParticles( &maxT, &numT );	
	printf("Total number of particles: %d\n", particlecount); 
	printf("Max Lifetime of a particle: %d\n", maxT ); 

	// Allocate memory for rendering quads
	if( triangleArray != NULL ) delete [] triangleArray;
	if( colorArray != NULL ) delete [] colorArray;
	if( normalArray != NULL ) delete [] normalArray;	

	triangleArray = new float[particlecount * 3 * 3];
	colorArray = new float[particlecount * 3 * 3];
	normalArray = new float[particlecount * 3 * 3];

	int binid = -1;
	float binlength = 0;
	if( sortEnabled == 1 )
	{
		// Get current modelview matrix
		_CurrentModelview();
	
		// Compute min-max z values in world space
		_ZMinMaxEyeSpace( z );
		printf("Z-range: [%f %f]\n", z[0], z[1]); 

		// Compute width of each bin
		binlength = (z[1] - z[0]) / (float)bincount;
	}

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
		int l = trace->size();

		int iP = 0;
		VECTOR3 p, nextp;
		for(list<VECTOR3*>::const_iterator
				pnIter = trace->begin(); 
			pnIter!= trace->end(); 
			pnIter++, iP++) 
		{
			// Get next particle
			p = **pnIter; 

			// If this is not the last point of the trace,
			// Get subsequent particle
			if(iP != (l-1))
			{
				list<VECTOR3*>::const_iterator tmpIter = pnIter;
				tmpIter++;
				nextp = **tmpIter;
			}
			else
			{
				nextp.Set( p(0), p(1), p(2) );
			}

			//printf( "P:%f, %f, %f\n", p(0), p(1), p(2) );
			//printf( "nextP:%f, %f, %f\n", nextp(0), nextp(1), nextp(2) );
						
			// Get World Space Coordinates of the particle
			VECTOR4 p4( p );
			p4 = modelview * p4; 

			if( sortEnabled == 1 )
			{
				// Compute bin id
				binid = (int)floor( p4(2) - z[0]) / binlength;
				if( binid<0 ) binid = 0;
				if( binid>=bincount ) binid = bincount - 1;
	
				// Compute array offset for this particle
				offset = (offsetToBins[binid] + pcountOfBins[binid]);

				// Shift array 
				int nParticlesToShift = 0;
				for(int i=binid+1; i<bincount; i++)
					nParticlesToShift += pcountOfBins[i];			
				memmove( (triangleArray+offset*9+9), (triangleArray+offset*9), nParticlesToShift * 9 * sizeof(float) );
				memmove( (colorArray+offset*9+9), (colorArray+offset*9), nParticlesToShift * 9 * sizeof(float) );
			}
	
			VECTOR3 dir = nextp - p;
			//dir.Normalize();
			
			float frac = 0.5f;
			float highx = p(0) + dir(0)*frac;//*g;
			float highy = p(1) + dir(1)*frac;//*g;
			float highz = p(2) + dir(2)*frac;//*g;
			
			//printf( "P:%f, %f, %f\n", p(0), p(1), p(2) );
			//printf( "nextP:%f, %f, %f\n", nextp(0), nextp(1), nextp(2) );
			//printf( "high:%f, %f, %f\n", highx, highy, highz );

			triangleArray[offset*9] = p(0);		triangleArray[offset*9+1] = p(1)+0.1f;	triangleArray[offset*9+2] = p(2);
			triangleArray[offset*9+3] = p(0);	triangleArray[offset*9+4] = p(1)-0.1f;	triangleArray[offset*9+5] = p(2);
			triangleArray[offset*9+6] = highx;	triangleArray[offset*9+7] = highy;		triangleArray[offset*9+8] = highz;

			// Calculate Normals
			if( lightEnabled == 1 )
			{
				VECTOR3 e1( triangleArray[offset*9+3]-triangleArray[offset*9],
					       triangleArray[offset*9+4]-triangleArray[offset*9+1],
						   triangleArray[offset*9+5]-triangleArray[offset*9+2] );
				VECTOR3 e2( triangleArray[offset*9+6]-triangleArray[offset*9+3],
					       triangleArray[offset*9+7]-triangleArray[offset*9+4],
						   triangleArray[offset*9+8]-triangleArray[offset*9+5] );		
				VECTOR3 e3( triangleArray[offset*9]-triangleArray[offset*9+6],
					       triangleArray[offset*9+1]-triangleArray[offset*9+7],
						   triangleArray[offset*9+2]-triangleArray[offset*9+8] );	
				
				e1.Normalize();
				e2.Normalize();
				e3.Normalize();

				VECTOR3 norm( e3(1)*e1(2) - e3(2)*e1(1), 
							  e3(2)*e1(0) - e3(0)*e1(2),
							  e3(0)*e1(1) - e3(1)*e1(0) );
				norm.Normalize();
				normalArray[offset*9] = norm(0);
				normalArray[offset*9+1] = norm(1);
				normalArray[offset*9+2] = norm(2);

				norm.Set( e1(1)*e2(2) - e1(2)*e2(1), 
							  e1(2)*e2(0) - e1(0)*e2(2),
							  e1(0)*e2(1) - e1(1)*e2(0) );
				norm.Normalize();
				normalArray[offset*9+3] = norm(0);
				normalArray[offset*9+4] = norm(1);
				normalArray[offset*9+5] = norm(2);

				norm.Set( e2(1)*e2(2) - e3(2)*e2(1), 
							  e2(2)*e2(0) - e3(0)*e2(2),
							  e2(0)*e2(1) - e3(1)*e2(0) );
				norm.Normalize();
				normalArray[offset*9+6] = norm(0);
				normalArray[offset*9+7] = norm(1);
				normalArray[offset*9+8] = norm(2);

				//VECTOR3 norm( (e1(0) + e2(0)) / 2.0f, (e1(1) + e2(1)) / 2.0f, (e1(2) + e2(2)) / 2.0f );
				//norm.Normalize();
		

			}

			if( sortEnabled == 1 )
			{
				// Update bin frequency
				pcountOfBins[binid] += 1;
			
				// Update bin offsets
				for(int i=binid+1; i<bincount; i++)
					offsetToBins[i] += 1;
			}

			if( currentColorScheme == 0 )
			{
				// Assign color based on time
				float f = (float)iP / (float)maxT;
				float f1 = (float)(iP+1) / (float)maxT;
				f = 0.2f + 0.8f * f;
				f1 = 0.2f + 0.8f * f1;
				colorArray[offset*9] = f;	colorArray[offset*9+1] = 0.4f; colorArray[offset*9+2] = 1-f;
				colorArray[offset*9+3] = f;	colorArray[offset*9+4] = 0.4f; colorArray[offset*9+5] = 1-f;
				colorArray[offset*9+6] = f1;	colorArray[offset*9+7] = 0.4f; colorArray[offset*9+8] = 1-f1;
			}
			else if( currentColorScheme == 1 )
			{
				// Assign color based on Glyph direction
				dir.Normalize();
				dir[0] = 0.2f + 0.8f * dir(0);
				dir[1] = 0.2f + 0.8f * dir(1);
				dir[2] = 0.2f + 0.8f * dir(2);

				colorArray[offset*9] = dir(0);		colorArray[offset*9+1] = dir(1); colorArray[offset*9+2] = dir(2);
				colorArray[offset*9+3] = dir(0);	colorArray[offset*9+4] = dir(1); colorArray[offset*9+5] = dir(2);
				colorArray[offset*9+6] = dir(0);	colorArray[offset*9+7] = dir(1); colorArray[offset*9+8] = dir(2);
			}

			// Increment offset, if no sorting
			if( sortEnabled == 0 )	offset ++;
		
		}// end inner for: for each trace

	}// end outer for

	delete [] z;
	delete [] offsetToBins;
	delete [] pcountOfBins;
}

void
CPointRendererGlyph::_TraversePoints_Glyph_gluSphere()
{
	printf( "Rendering Spherical Glyphs with lighting ...\n" );
	
	// Initialize local variables
	const list<vtListSeedTrace*>* sl_list = (const list<vtListSeedTrace*>*)this->pDataSource;
	int iT=0;
	int bincount = 20;
	float *z = new float[2];
	int *offsetToBins = new int[bincount];
	int *pcountOfBins = new int[bincount];
	for(int i=0; i<bincount; i++)
	{
		offsetToBins[i] = 0;
		pcountOfBins[i] = 0;
	}

	// Get total number of particles
	// Also, maximum lifetime of a particle
	maxT = 0;
	particlecount = _CountParticles( &maxT, &numT );	
	printf("Total number of particles: %d\n", particlecount); 
	printf("Max Lifetime of a particle: %d\n", maxT ); 

	// Allocate memory for rendering quads
	if( pointArray != NULL ) delete [] pointArray;
	if( colorArray != NULL ) delete [] colorArray;

	pointArray = new float[particlecount * 3];
	colorArray = new float[particlecount * 3];

	int binid = -1;
	float binlength = 0;
	if( sortEnabled == 1 )
	{
		// Get current modelview matrix
		_CurrentModelview();
	
		// Compute min-max z values in world space
		_ZMinMaxEyeSpace( z );
		printf("Z-range: [%f %f]\n", z[0], z[1]); 

		// Compute width of each bin
		binlength = (z[1] - z[0]) / (float)bincount;
	}

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

			if( sortEnabled == 1 )
			{
				// Compute bin id
				binid = (int)floor( p4(2) - z[0]) / binlength;
				if( binid<0 ) binid = 0;
				if( binid>=bincount ) binid = bincount - 1;
					
				// Compute array offset for this particle
				offset = (offsetToBins[binid] + pcountOfBins[binid]);

				// Shift array 
				int nParticlesToShift = 0;
				for(int i=binid+1; i<bincount; i++)
					nParticlesToShift += pcountOfBins[i];			
				memmove( (pointArray+offset*3+3), (pointArray+offset*3), nParticlesToShift * 3 * sizeof(float) );
				memmove( (colorArray+offset*3+3), (colorArray+offset*3), nParticlesToShift * 3 * sizeof(float) );
			}
			
			// Store vertex location
			pointArray[offset*3] = p(0);	pointArray[offset*3+1] = p(1); pointArray[offset*3+2] = p(2);			
			
			if( currentColorScheme == 0 )
			{
				// Assign color based on time
				float f = (float)iP / (float)maxT;
				colorArray[offset*3] = f;	colorArray[offset*3+1] = 0.5f; colorArray[offset*3+2] = 1-f;
			}


			if( sortEnabled == 1)
			{
				// Update Count
				pcountOfBins[binid] += 1;
			
				// Update offsets
				for(int i=binid+1; i<bincount; i++)
					offsetToBins[i] += 1;
			}
			else
				offset ++;
		
		}// end inner for: for each trace

	}// end outer for

	delete [] z;
	delete [] offsetToBins;
	delete [] pcountOfBins;

}// end function

void
CPointRendererGlyph::_TraversePoints_Glyph_gluCone()
{
	printf( "Rendering Cone-shaped Glyphs with lighting...\n" );
	
	// Initialize local variables
	const list<vtListSeedTrace*>* sl_list = (const list<vtListSeedTrace*>*)this->pDataSource;
	int iT=0;
	int bincount = 20;
	float *z = new float[2];
	int *offsetToBins = new int[bincount];
	int *pcountOfBins = new int[bincount];
	for(int i=0; i<bincount; i++)
	{
		offsetToBins[i] = 0;
		pcountOfBins[i] = 0;
	}

	// Get total number of particles to begin
	maxT = 0;
	particlecount = _CountParticles( &maxT, &numT );	
	printf("Total number of particles: %d\n", particlecount); 
	printf("Max Lifetime of a particle: %d\n", maxT ); 

	// Allocate memory for rendering cylinders
	if( pointArray != NULL ) delete [] pointArray;
	if( colorArray != NULL ) delete [] colorArray;
	if( normalArray != NULL ) delete [] normalArray;
	if( angleArray != NULL ) delete [] angleArray;
	if( heightArray != NULL ) delete [] heightArray;

	pointArray = new float[particlecount * 3];
	colorArray = new float[particlecount * 3];
	normalArray = new float[particlecount * 3];
	angleArray = new float[particlecount];
	heightArray = new float[particlecount];

	int binid = -1;
	float binlength = 0;
	if( sortEnabled == 1 )
	{
		// Get current modelview matrix
		_CurrentModelview();
	
		// Compute min-max z values in world space
		_ZMinMaxEyeSpace( z );
		printf("Z-range: [%f %f]\n", z[0], z[1]); 

		// Compute width of each bin
		binlength = (z[1] - z[0]) / (float)bincount;
	}

	
	iT = 0;
	int offset = 0;	
	float frac = 0.1f;
	VECTOR3 zaxis( 0, 0, 1 );

	// Scan through all points and asssign bins to them
	for(list<vtListSeedTrace*>::const_iterator
			pIter =  sl_list->begin(); 
		pIter!=sl_list->end(); 
		pIter++, iT++) 
	{
	    // Get next trace
		const vtListSeedTrace *trace = *pIter; 
		int l = trace->size();

		int iP = 0;
		VECTOR3 p, nextp;
		for(list<VECTOR3*>::const_iterator
				pnIter = trace->begin(); 
			pnIter!= trace->end(); 
			pnIter++, iP++) 
		{
			// Get next particle
			p = **pnIter; 

			// If this is not the last point of the trace,
			// Get subsequent particle
			if(iP != (l-1))
			{
				list<VECTOR3*>::const_iterator tmpIter = pnIter;
				tmpIter++;
				nextp = **tmpIter;
			}
			else
			{
				nextp.Set( p(0), p(1), p(2) );
			}

			//printf( "P:%f, %f, %f\n", p(0), p(1), p(2) );
			//printf( "nextP:%f, %f, %f\n", nextp(0), nextp(1), nextp(2) );
						
			// Get World Space Coordinates of the particle
			VECTOR4 p4( p );
			p4 = modelview * p4; 

			if( sortEnabled == 1 )
			{
				// Compute bin id
				binid = (int)floor( p4(2) - z[0]) / binlength;
				if( binid<0 ) binid = 0;
				if( binid>=bincount ) binid = bincount - 1;
	
				// Compute array offset for this particle
				offset = (offsetToBins[binid] + pcountOfBins[binid]);

				// Shift array 
				int nParticlesToShift = 0;
				for(int i=binid+1; i<bincount; i++)
					nParticlesToShift += pcountOfBins[i];			
				memmove( (triangleArray+offset*9+9), (triangleArray+offset*9), nParticlesToShift * 9 * sizeof(float) );
				memmove( (colorArray+offset*9+9), (colorArray+offset*9), nParticlesToShift * 9 * sizeof(float) );
			}
	
			VECTOR3 dir = nextp - p;
							
			float highx = p(0) + dir(0)*frac;
			float highy = p(1) + dir(1)*frac;
			float highz = p(2) + dir(2)*frac;
			
			dir.Normalize();
			
			VECTOR3 norm( zaxis(1)*dir(2) - zaxis(2)*dir(1), 
						  zaxis(2)*dir(0) - zaxis(0)*dir(2),
						  zaxis(0)*dir(1) - zaxis(1)*dir(0) );
			//norm.Normalize();

			//printf( "P:%f, %f, %f\n", p(0), p(1), p(2) );
			//printf( "nextP:%f, %f, %f\n", nextp(0), nextp(1), nextp(2) );
			//printf( "high:%f, %f, %f\n", highx, highy, highz );

			pointArray[offset*3] = p(0);	pointArray[offset*3+1] = p(1);	pointArray[offset*3+2] = p(2);
			normalArray[offset*3] = norm(0);	normalArray[offset*3+1] = norm(1);	normalArray[offset*3+2] = norm(2);
			angleArray[offset] =  ( 180.0f / 3.14f ) * acos( dir(2) ); 
			heightArray[offset] = sqrt( ( highx - p(0) )*(highx - p(0) ) + 
										( highy - p(1) )*(highy - p(1) ) + 
										( highy - p(2) )*(highy - p(2) ) );

			if( sortEnabled == 1 )
			{
				// Update bin frequency
				pcountOfBins[binid] += 1;
			
				// Update bin offsets
				for(int i=binid+1; i<bincount; i++)
					offsetToBins[i] += 1;
			}

			if( currentColorScheme == 0 )
			{
				// Assign color based on time
				float f = (float)iP / (float)maxT;
				float f1 = (float)(iP+1) / (float)maxT;
				f = 0.2f + 0.8f * f;
				f1 = 0.2f + 0.8f * f1;
				colorArray[offset*3] = f;	colorArray[offset*3+1] = 0.4f; colorArray[offset*3+2] = 1-f;
			}
			else if( currentColorScheme == 1 )
			{
				// Assign color based on Glyph direction
				dir.Normalize();
				dir[0] = 0.2f + 0.8f * dir(0);
				dir[1] = 0.2f + 0.8f * dir(1);
				dir[2] = 0.2f + 0.8f * dir(2);

				colorArray[offset*3] = dir(0);		colorArray[offset*3+1] = dir(1); colorArray[offset*3+2] = dir(2);
			}

			// Increment offset, if no sorting
			if( sortEnabled == 0 )	offset ++;
		
		}// end inner for: for each trace

	}// end outer for


	delete [] z;
	delete [] offsetToBins;
	delete [] pcountOfBins;
}

void
CPointRendererGlyph::_RenderPoints_Glyph_Arrow( )
{
	if( lightEnabled == 1)
	{
		glEnable( GL_LIGHTING );
		glEnable( GL_LIGHT0 );
		glEnable( GL_DEPTH_TEST );
		
		//glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT );
		//glEnable( GL_COLOR_MATERIAL );
		//glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );

		GLfloat mat_ambient[] = { 0.1, 0.1, 0.1, 1 };
		GLfloat mat_diffuse[] = { 0.5, 0.5, 0.5, 1 };
		GLfloat mat_specular[] = { 1, 1, 1, 1 };
		GLfloat mat_emission[] = { 0.2, 0.2, 0.2, 1 };
		GLfloat mat_shininess[] = { 50 };

	    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient );
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission); 
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	}

	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	glEnable (GL_LINE_SMOOTH);

	glClearColor (0.0, 0.0, 0.0, 0.0);
	glShadeModel (GL_SMOOTH);

	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_COLOR_ARRAY );
	if( lightEnabled == 1 )	glEnableClientState( GL_NORMAL_ARRAY );
	
	glVertexPointer( 3, GL_FLOAT, 0, triangleArray );
	glColorPointer( 3, GL_FLOAT, 0, colorArray );
	if( lightEnabled == 1 )	glNormalPointer( GL_FLOAT, 0, normalArray );
	
	// Render the vertex array	
	glDrawArrays( GL_TRIANGLES, 0, particlecount * 3 ) ;
	//glutSolidTeapot(10.0f);

	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_COLOR_ARRAY );
	
	if( lightEnabled == 1)
	{
		glDisableClientState( GL_NORMAL_ARRAY );
		glDisable( GL_LIGHTING );
		glDisable( GL_LIGHT0 );
	}

}// end function

void 
CPointRendererGlyph::_RenderPoints_Glyph_gluSphere()
{
	if( lightEnabled == 1 )
	{
		glEnable( GL_LIGHTING );
		glEnable( GL_LIGHT0 );
		glEnable( GL_DEPTH_TEST );
		
		//glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT );
		//glEnable( GL_COLOR_MATERIAL );
	
		GLfloat mat_ambient[] = { 0.1, 0.1, 0.1, 1 };
		GLfloat mat_diffuse[] = { 0.3, 0.3, 0.3, 1 };
		GLfloat mat_specular[] = { 1, 1, 1, 1 };
		GLfloat mat_emission[] = { 0, 0, 0, 1 };
		GLfloat mat_shininess[] = { 50 };

	    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient );
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission); 
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
	}

	//glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	//glEnable (GL_LINE_SMOOTH);

	glClearColor (0.0, 0.0, 0.0, 0.0);
	glShadeModel (GL_SMOOTH);

	GLUquadric *myQuad = gluNewQuadric();
	GLdouble radius = 0.2f;
	GLint slices=20, stacks=10;

	gluQuadricNormals( myQuad, GLU_SMOOTH );
	gluQuadricOrientation( myQuad, GLU_OUTSIDE );
			
	for( int i=0; i<particlecount; i++ )
	{
		glPushMatrix();
		glTranslatef( pointArray[i*3], pointArray[i*3+1], pointArray[i*3+2] );
		gluSphere( myQuad , radius , slices , stacks  );
		glColor3f( colorArray[i*3], colorArray[i*3+1], colorArray[i*3+2] );
		glPopMatrix();
	}// end

	gluDeleteQuadric( myQuad );	

}// end function

void 
CPointRendererGlyph::_RenderPoints_Glyph_gluCone()
{
	if( lightEnabled == 1 )
	{
		//glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT );
		//glEnable( GL_COLOR_MATERIAL );
		//glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );

		GLfloat mat_ambient[] = { 0.1, 0.1, 0.1, 1 };
		GLfloat mat_diffuse[] = { 0.1, 0.1, 0.1, 1 };
		GLfloat mat_specular[] = { 0.1, 0.1, 0.1, 1 };
		GLfloat mat_emission[] = { 0, 0, 0, 1 };
		GLfloat mat_shininess[] = { 50 };

	    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient );
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission); 
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

		glEnable( GL_LIGHTING );
		glEnable( GL_LIGHT0 );
		glEnable( GL_DEPTH_TEST );
	}

	//glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
	//glEnable (GL_LINE_SMOOTH);

	glClearColor (0.0, 0.0, 0.0, 0.0);
	glShadeModel (GL_SMOOTH);

	GLUquadric *myQuad = gluNewQuadric();
	GLdouble base = 0.1f;
	GLdouble top = 0.05;
	GLdouble height = 1.5f;
	GLint slices=30, stacks=10;

	gluQuadricNormals( myQuad, GLU_SMOOTH );
	gluQuadricOrientation( myQuad, GLU_OUTSIDE );

	for( int i=0; i<particlecount; i++ )
	{
		glPushMatrix();
		glRotatef( angleArray[i], normalArray[i*3], normalArray[i*3+1], normalArray[i*3+2] );
		glTranslatef( pointArray[i*3], pointArray[i*3+1], pointArray[i*3+2] );
		//gluCylinder( myQuad , base, top, heightArray[i], slices, stacks  );
		gluCylinder( myQuad , base, top, height, slices, stacks  );
		//glColor3f( colorArray[i*3], colorArray[i*3+1], colorArray[i*3+2] );
		glPopMatrix();
	}// end

	if( lightEnabled == 1)
	{
		glDisable( GL_LIGHTING );
		glDisable( GL_LIGHT0 );
	}

	gluDeleteQuadric( myQuad );	

}// end function

void CPointRendererGlyph::_Draw()
{
	if( 0 == uPid )
		return;

	glCallList(uPid);

}// end function

void CPointRendererGlyph::_Update()
{
	if( NULL == pDataSource )
		return;

	if( 0 == uPid )
		uPid = glGenLists(1);

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	if( currentGlyphType == GLYPH_ARROW )
	{
		_TraversePoints_Glyph_Arrow();
		_RenderPoints_Glyph_Arrow();		
	}
	else if( currentGlyphType == GLYPH_SPHERE )
	{
		_TraversePoints_Glyph_gluSphere();
		_RenderPoints_Glyph_gluSphere();
	}
	else if( currentGlyphType == GLYPH_CONE )
	{
		_TraversePoints_Glyph_gluCone();
		_RenderPoints_Glyph_gluCone();
	}

	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();

}// end function

void CPointRendererGlyph::_UpdateGlyphType( int inGlyphType  )
{
	if( NULL == pDataSource )
		return;

	if( currentGlyphType != inGlyphType )
	{
		currentGlyphType = inGlyphType;
		uPid = glGenLists(1);
	}
	else
		return;

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	if( currentGlyphType == GLYPH_ARROW )
	{
		_TraversePoints_Glyph_Arrow();
		_RenderPoints_Glyph_Arrow();
	}
	else if( currentGlyphType == GLYPH_SPHERE )
	{
		_TraversePoints_Glyph_gluSphere();
		_RenderPoints_Glyph_gluSphere();
	}
	else if( currentGlyphType == GLYPH_CONE )
	{
		_TraversePoints_Glyph_gluCone();
		_RenderPoints_Glyph_gluCone();
	}
	
	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();

}// end function

void CPointRendererGlyph::_UpdateSorting( int sortFlag )
{
	if( NULL == pDataSource )
		return;

	if( sortEnabled != sortFlag )
	{
		sortEnabled = sortFlag;
		uPid = glGenLists(1);
	}
	else
		return;

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	if( currentGlyphType == GLYPH_ARROW )		_TraversePoints_Glyph_Arrow();
	else if( currentGlyphType == GLYPH_SPHERE )	_TraversePoints_Glyph_gluSphere();
	else if( currentGlyphType == GLYPH_CONE )	_TraversePoints_Glyph_gluCone();

	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();

}// end function

void CPointRendererGlyph::_UpdateColorScheme( int inColorScheme )
{
	if( NULL == pDataSource )
		return;

	if( currentColorScheme != inColorScheme )
	{
		currentColorScheme = inColorScheme;
		uPid = glGenLists(1);
	}
	else
		return;

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	if( currentGlyphType == GLYPH_ARROW )		_RenderPoints_Glyph_Arrow();
	else if( currentGlyphType == GLYPH_SPHERE )	_RenderPoints_Glyph_gluSphere();
	else if( currentGlyphType == GLYPH_CONE )	_RenderPoints_Glyph_gluCone();

	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();

}// end function

void CPointRendererGlyph::_UpdateLighting( int lightFlag )
{
	if( NULL == pDataSource )
		return;

	if( lightEnabled != lightFlag )
	{
		lightEnabled = lightFlag;
		uPid = glGenLists(1);
	}
	else
		return;

	glNewList(uPid, GL_COMPILE);

	// Apply transformations
	_TraversePointsBegin();

	if( currentGlyphType == GLYPH_ARROW )		_RenderPoints_Glyph_Arrow();
	else if( currentGlyphType == GLYPH_SPHERE )	_RenderPoints_Glyph_gluSphere();
	else if( currentGlyphType == GLYPH_CONE )	_RenderPoints_Glyph_gluCone();

	// Unapply Transformation
	_TraversePointsEnd();

	glEndList();

}// end function

CPointRendererGlyph::CPointRendererGlyph(void)
{
	// Initialize data structure
	pointArray = NULL;
	colorArray = NULL;
	normalArray = NULL;
	triangleArray = NULL;
	angleArray = NULL;
	heightArray = NULL;

	currentGlyphType = GLYPH_SPHERE;
	sortEnabled = 0;
	currentColorScheme = 0;
	lightEnabled = 1;


}// end function

CPointRendererGlyph::~CPointRendererGlyph(void)
{
}// end function
