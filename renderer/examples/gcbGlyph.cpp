/*

gcbGlyphs:
The is a demo to show how to render the particles generated from flow field as small 3D glyphs.

*/

#include <list>
#include <iterator>

#include "gcb.h"
#include "OSUFlow.h"
#include "PointRendererGlyph.h"

OSUFlow *osuflow; 
VECTOR3 minLen, maxLen; 
list<vtListSeedTrace*> sl_list; 
float center[3], len[3]; 

CPointRendererGlyph cPointRendererGlyph;

////////////////////////////////////////////////////////////////////////////
void compute_streamlines() 
{
  LOG("");

  float from[3], to[3]; 

  from[0] = minLen[0];   from[1] = minLen[1];   from[2] = minLen[2]; 
  to[0] = maxLen[0];   to[1] = maxLen[1];   to[2] = maxLen[2]; 

  printf("generating seeds...\n"); 
  osuflow->SetRandomSeedPoints(from, to, 100); 
  int nSeeds; 
  VECTOR3* seeds = osuflow->GetSeeds(nSeeds); 
  for (int i=0; i<nSeeds; i++) 
    printf(" seed no. %d : [%f %f %f]\n", i, seeds[i][0], seeds[i][1], seeds[i][2]); 

  sl_list.clear(); 

  printf("compute streamlines..\n"); 
  osuflow->SetIntegrationParams(1, 5); 
  osuflow->GenStreamLines(sl_list , FORWARD_DIR, 500, 0); 
  printf(" done integrations\n"); 
  printf("list size = %d\n", (int)sl_list.size()); 

  cPointRendererGlyph._Update();
}

void draw_glyphs() 
{
	cPointRendererGlyph._Draw();
}

///////////////////////////////////////////////////////////////////////////////
void
_KeyboardFunc(unsigned char ubKey, int iX, int iY)
{
	switch(ubKey)
	{
	case 's':
		compute_streamlines();
		glutPostRedisplay();
		break;

	case 'q':
		cPointRendererGlyph._UpdateGlyphType( 1 );
		glutPostRedisplay();
		break;
	
	case 'w':
		cPointRendererGlyph._UpdateGlyphType( 2 );
		glutPostRedisplay();
		break;
	
	case 'e':
		cPointRendererGlyph._UpdateGlyphType( 3 );
		glutPostRedisplay();
		break;
	
	}

}

void
_DisplayFunc()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// render the scene
    draw_glyphs(); 

	// NOTE: Call glutSwapBuffers() at the end of your display function
	glutSwapBuffers();
}



void
init()
{
	LOG(printf("Initialize here."));
	glEnable(GL_DEPTH_TEST);
}

void 
quit()
{
	LOG(printf("Clean up here."));
}

int
main(int argn, char* argv[])
{
	///////////////////////////////////////////////////////////////
	// when use GCB, it is still needed to initialize GLUT
	glutInit(&argn, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_STENCIL );

	///////////////////////////////////////////////////////////////
	// initialize OSU flow
	osuflow = new OSUFlow(); 

	// load the scalar field
	LOG(printf("read file %s\n", argv[1]));
	osuflow->LoadData((const char*)argv[1], true); //true: a steady flow field 
	// comptue the bounding box of the streamlines 
	VECTOR3 minB, maxB; 
	osuflow->Boundary(minLen, maxLen); // get the boundary 
	minB[0] = minLen[0]; minB[1] = minLen[1];  minB[2] = minLen[2];
	maxB[0] = maxLen[0]; maxB[1] = maxLen[1];  maxB[2] = maxLen[2];
	//  osuflow->SetBoundary(minB, maxB);  // set the boundary. just to test
										 // the subsetting feature of OSUFlow
	printf("volume boundary X: [%f %f] Y: [%f %f] Z: [%f %f]\n", 
								minLen[0], maxLen[0], minLen[1], maxLen[1], 
								minLen[2], maxLen[2]); 

	center[0] = (minLen[0]+maxLen[0])/2.0; 
	center[1] = (minLen[1]+maxLen[1])/2.0; 
	center[2] = (minLen[2]+maxLen[2])/2.0; 
	printf("center is at %f %f %f \n", center[0], center[1], center[2]); 
	len[0] = maxLen[0]-minLen[0]; 
	len[1] = maxLen[1]-minLen[1]; 
	len[2] = maxLen[2]-minLen[2]; 

	///////////////////////////////////////////////////////////////
	cPointRendererGlyph._SetBoundingBox(
		minLen[0], minLen[1], minLen[2], 
		maxLen[0], maxLen[1], maxLen[2]);
	cPointRendererGlyph._SetDataSource(&sl_list);
	
	///////////////////////////////////////////////////////////////
	glutCreateWindow("GCB Glyphs");

	// specify the callbacks you really need. Except gcbInit() and gcbDisplayFunc(), other callbacks are optional
	gcbInit(init, quit);
	gcbDisplayFunc(_DisplayFunc);
	gcbKeyboardFunc(_KeyboardFunc);

	// enter the GLUT loop
	glutMainLoop();

	return 0;
}

/*

$Log: gcbGlyphs.cpp,v $
Revision 1.1  2010/04/19 16:41:38  chaudhua
*** empty log message ***


*/
