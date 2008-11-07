//////////////////////////////////////////////////////////////
//
//  Triangulation functions shared by several isosurface classes
//
//
//////////////////////////////////////////////////////////////

#include "mcube.h"
#include "triangulator.h"

//#define gl_vertex(f) glVertex3f(f[0], f[1], f[2]); 
//#define gl_normal(n) glNormal3f(n[0], n[1], n[2]); 

int vt_off[9]; 

///////////////////////////////////////////////////////
//
//   Compute the normal -- just a simple cross-product. 
//   If smoother surfaces are required, better calculate
//   the gradient from the data field. 
//
inline int compute_normal(VECTOR3 f1, VECTOR3 f2, VECTOR3 f3, 
						  VECTOR3 &normal) 
{
	float x1 = f2[0] - f1[0]; 
	float y1 = f2[1] - f1[1]; 
	float z1 = f2[2] - f1[2]; 

	float x2 = f3[0] - f2[0]; 
	float y2 = f3[1] - f2[1]; 
	float z2 = f3[2] - f2[2]; 

	normal[0] = y1*z2 - y2*z1; 
	normal[1] = x2*z1 - x1*z2; 
	normal[2] = x1*y2 - x2*y1; 

	float mag = (float) sqrt( (double)(normal[0] * normal[0] + 
		normal[1] * normal[1] + 
		normal[2] * normal[2])); 

	if( fabs(mag) < 1.0e-06)
		return -1;

	normal[0] = normal[0]/mag; 
	normal[1] = normal[1]/mag; 
	normal[2] = normal[2]/mag; 
	return 1;
}

////////////////////////////////////////////////////////

int iso_hex(int index, int xdim, int ydim, 
			int zdim, float isoval, float* field,
			vector<VECTOR3*>& vTriangles)
{
//	glBegin(GL_TRIANGLES); 
	float oval[9];
	float ov[9][3];     

	int i = index % (xdim-1); 
	int j = (index / (xdim-1)) % (ydim-1); 
	int k = index / ((xdim-1)*(ydim-1)); 

	int idx = k*xdim*ydim+j*xdim+i; 

	oval[1]=field[idx]-isoval;
	oval[2]=field[idx+vt_off[2]]-isoval;
	oval[3]=field[idx+vt_off[3]]-isoval;
	oval[4]=field[idx+vt_off[4]]-isoval;
	oval[5]=field[idx+vt_off[5]]-isoval;
	oval[6]=field[idx+vt_off[6]]-isoval;
	oval[7]=field[idx+vt_off[7]]-isoval;
	oval[8]=field[idx+vt_off[8]]-isoval;

	ov[1][0] = i; ov[1][1] = j; ov[1][2] = k; 
	ov[2][0] = i+1; ov[2][1] = j; ov[2][2] = k; 
	ov[3][0] = i+1; ov[3][1] = j+1; ov[3][2] = k; 
	ov[4][0] = i; ov[4][1] = j+1; ov[4][2] = k; 
	ov[5][0] = i; ov[5][1] = j; ov[5][2] = k+1; 
	ov[6][0] = i+1; ov[6][1] = j; ov[6][2] = k+1; 
	ov[7][0] = i+1; ov[7][1] = j+1; ov[7][2] = k+1; 
	ov[8][0] = i; ov[8][1] = j+1; ov[8][2] = k+1; 

	int mask=0;
	for(int index1=1; index1<=8; index1++){
		if(oval[index1]<0)
			mask|=1<<(index1-1);
	}
	MCubeTable* tab=&mcube_table[mask];
	double val[9];
	float  v[9][3]; 
	for(int indx=1;indx<=8;indx++){
		val[indx]=oval[tab->permute[indx-1]];
		for (int i=0; i<3; i++) 
			v[indx][i]=ov[tab->permute[indx-1]][i];
	}
	int wcase=tab->which_case;
	int no_faces = 0; 
	VECTOR3 f1,f2, f3, f4, f5, f6; 
	VECTOR3 f7,f8, f9, f10, f11, f12; 
	VECTOR3 normal; 


	switch(wcase){
	case 0:
		break;

	case 1:
		{
			float weight = val[1]/(val[1]-val[2]); 

			f1[0] = v[1][0] * (1-weight) + v[2][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[2][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[2][2] * weight; 

			weight = val[1]/(val[1]-val[5]); 
			f2[0] = v[1][0] * (1-weight) + v[5][0] * weight; 
			f2[1] = v[1][1] * (1-weight) + v[5][1] * weight; 
			f2[2] = v[1][2] * (1-weight) + v[5][2] * weight; 

			weight = val[1]/(val[1]-val[4]); 
			f3[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f3[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f3[2] = v[1][2] * (1-weight) + v[4][2] * weight;

			compute_normal(f1, f2, f3, normal); 

			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));
						
			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1); 
			geom->vertex3fv(f2); 
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/
		}  
		no_faces = 1; 
		break;
	case 2:
		{
			double weight = val[1]/(val[1]-val[5]); 
			f1[0] = v[1][0] * (1-weight) + v[5][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[5][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[5][2] * weight; 

			weight = val[2]/(val[2]-val[6]); 
			f2[0] = v[2][0] * (1-weight) + v[6][0] * weight; 
			f2[1] = v[2][1] * (1-weight) + v[6][1] * weight; 
			f2[2] = v[2][2] * (1-weight) + v[6][2] * weight; 

			weight = val[2]/(val[2]-val[3]); 
			f3[0] = v[2][0] * (1-weight) + v[3][0] * weight; 
			f3[1] = v[2][1] * (1-weight) + v[3][1] * weight; 
			f3[2] = v[2][2] * (1-weight) + v[3][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);  geom->vertex3fv(f2);  geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[1]/(val[1]-val[4]); 
			f4[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f4[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f4[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f3, f4, f1, normal); 
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f3, f4, f1, normal) != -1 ){
			geom->vertex3fv(f3);  geom->vertex3fv(f4);  geom->vertex3fv(f1); 
			geom->normal3fv( normal); 
			// 3,4,1
			}
			*/
		}
		no_faces = 2; 
		break;
	case 3:
		{
			double weight = val[1]/(val[1]-val[2]); 
			f1[0] = v[1][0] * (1-weight) + v[2][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[2][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[2][2] * weight; 

			weight = val[1]/(val[1]-val[5]); 
			f2[0] = v[1][0] * (1-weight) + v[5][0] * weight; 
			f2[1] = v[1][1] * (1-weight) + v[5][1] * weight; 
			f2[2] = v[1][2] * (1-weight) + v[5][2] * weight; 

			weight = val[1]/(val[1]-val[4]); 
			f3[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f3[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f3[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);   geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[3]/(val[3]-val[2]); 
			f4[0] = v[3][0] * (1-weight) + v[2][0] * weight; 
			f4[1] = v[3][1] * (1-weight) + v[2][1] * weight; 
			f4[2] = v[3][2] * (1-weight) + v[2][2] * weight; 

			weight = val[3]/(val[3]-val[7]); 
			f5[0] = v[3][0] * (1-weight) + v[7][0] * weight; 
			f5[1] = v[3][1] * (1-weight) + v[7][1] * weight; 
			f5[2] = v[3][2] * (1-weight) + v[7][2] * weight; 

			weight = val[3]/(val[3]-val[4]); 
			f6[0] = v[3][0] * (1-weight) + v[4][0] * weight; 
			f6[1] = v[3][1] * (1-weight) + v[4][1] * weight; 
			f6[2] = v[3][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f4, f5, f6, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f4, f5, f6, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f5);   geom->vertex3fv(f6); 
			geom->normal3fv( normal); 
			// 4,5,6
			}
			*/
		}
		no_faces = 2; 
		break;
	case 4:
		{
			double weight = val[1]/(val[1]-val[2]); 
			f1[0] = v[1][0] * (1-weight) + v[2][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[2][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[2][2] * weight; 

			weight = val[1]/(val[1]-val[5]); 
			f2[0] = v[1][0] * (1-weight) + v[5][0] * weight; 
			f2[1] = v[1][1] * (1-weight) + v[5][1] * weight; 
			f2[2] = v[1][2] * (1-weight) + v[5][2] * weight; 

			weight = val[1]/(val[1]-val[4]); 
			f3[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f3[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f3[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);   geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[7]/(val[7]-val[3]); 

			f4[0] = v[7][0] * (1-weight) + v[3][0] * weight; 
			f4[1] = v[7][1] * (1-weight) + v[3][1] * weight; 
			f4[2] = v[7][2] * (1-weight) + v[3][2] * weight; 

			weight = val[7]/(val[7]-val[8]); 
			f5[0] = v[7][0] * (1-weight) + v[8][0] * weight; 
			f5[1] = v[7][1] * (1-weight) + v[8][1] * weight; 
			f5[2] = v[7][2] * (1-weight) + v[8][2] * weight; 

			weight = val[7]/(val[7]-val[6]); 
			f6[0] = v[7][0] * (1-weight) + v[6][0] * weight; 
			f6[1] = v[7][1] * (1-weight) + v[6][1] * weight; 
			f6[2] = v[7][2] * (1-weight) + v[6][2] * weight; 

			compute_normal(f4, f5, f6, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f4, f5, f6, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f5);   geom->vertex3fv(f6); 
			geom->normal3fv( normal); 
			// 4,5,6
			}
			*/
		}
		no_faces = 2; 
		break;
	case 5:
		{
			double weight = val[2]/(val[2]-val[1]); 
			f1[0] = v[2][0] * (1-weight) + v[1][0] * weight; 
			f1[1] = v[2][1] * (1-weight) + v[1][1] * weight; 
			f1[2] = v[2][2] * (1-weight) + v[1][2] * weight; 

			weight = val[2]/(val[2]-val[3]); 
			f2[0] = v[2][0] * (1-weight) + v[3][0] * weight; 
			f2[1] = v[2][1] * (1-weight) + v[3][1] * weight; 
			f2[2] = v[2][2] * (1-weight) + v[3][2] * weight; 

			weight = val[5]/(val[5]-val[1]); 
			f3[0] = v[5][0] * (1-weight) + v[1][0] * weight; 
			f3[1] = v[5][1] * (1-weight) + v[1][1] * weight; 
			f3[2] = v[5][2] * (1-weight) + v[1][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);   geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[5]/(val[5]-val[8]); 
			f4[0] = v[5][0] * (1-weight) + v[8][0] * weight; 
			f4[1] = v[5][1] * (1-weight) + v[8][1] * weight; 
			f4[2] = v[5][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f4, f3, f2, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f4, f3, f2, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f3);   geom->vertex3fv(f2); 
			geom->normal3fv( normal); 
			// 4,3,2
			}
			*/
			weight = val[6]/(val[6]-val[7]); 
			f5[0] = v[6][0] * (1-weight) + v[7][0] * weight; 
			f5[1] = v[6][1] * (1-weight) + v[7][1] * weight; 
			f5[2] = v[6][2] * (1-weight) + v[7][2] * weight; 

			compute_normal(f5, f4, f2, normal); 
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f5, f4, f2, normal) != -1 ){
			geom->vertex3fv(f5);
			geom->vertex3fv(f4);   geom->vertex3fv(f2); 
			geom->normal3fv( normal); 
			// 5,4,2
			}
			*/
		}
		no_faces = 3; 
		break;
	case 6:
		{
			double weight = val[1]/(val[1]-val[5]); 
			f1[0] = v[1][0] * (1-weight) + v[5][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[5][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[5][2] * weight; 

			weight = val[2]/(val[2]-val[6]); 
			f2[0] = v[2][0] * (1-weight) + v[6][0] * weight; 
			f2[1] = v[2][1] * (1-weight) + v[6][1] * weight; 
			f2[2] = v[2][2] * (1-weight) + v[6][2] * weight; 

			weight = val[2]/(val[2]-val[3]); 
			f3[0] = v[2][0] * (1-weight) + v[3][0] * weight; 
			f3[1] = v[2][1] * (1-weight) + v[3][1] * weight; 
			f3[2] = v[2][2] * (1-weight) + v[3][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);   geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[1]/(val[1]-val[4]); 

			f4[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f4[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f4[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f3, f4, f1, normal); 
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(normal));


			/*
			if( compute_normal(f3, f4, f1, normal) != -1 ){
			geom->vertex3fv(f3);
			geom->vertex3fv(f4);   geom->vertex3fv(f1); 
			geom->normal3fv( normal); 
			// 3,4,1
			}
			*/

			weight = val[7]/(val[7]-val[3]); 
			f5[0] = v[7][0] * (1-weight) + v[3][0] * weight; 
			f5[1] = v[7][1] * (1-weight) + v[3][1] * weight; 
			f5[2] = v[7][2] * (1-weight) + v[3][2] * weight; 

			weight = val[7]/(val[7]-val[8]); 
			f6[0] = v[7][0] * (1-weight) + v[8][0] * weight; 
			f6[1] = v[7][1] * (1-weight) + v[8][1] * weight; 
			f6[2] = v[7][2] * (1-weight) + v[8][2] * weight; 

			weight = val[7]/(val[7]-val[6]); 
			f7[0] = v[7][0] * (1-weight) + v[6][0] * weight; 
			f7[1] = v[7][1] * (1-weight) + v[6][1] * weight; 
			f7[2] = v[7][2] * (1-weight) + v[6][2] * weight; 

			compute_normal(f5, f6, f7, normal); 
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(f7));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f5, f6, f7, normal) != -1 ){
			geom->vertex3fv(f5);
			geom->vertex3fv(f6);   geom->vertex3fv(f7); 
			geom->normal3fv( normal); 
			// 5,6,7
			}
			*/
		}
		no_faces = 3; 
		break;
	case 7:
		{
			double weight = val[2]/(val[2]-val[1]); 
			f1[0] = v[2][0] * (1-weight) + v[1][0] * weight; 
			f1[1] = v[2][1] * (1-weight) + v[1][1] * weight; 
			f1[2] = v[2][2] * (1-weight) + v[1][2] * weight; 

			weight = val[2]/(val[2]-val[3]); 
			f2[0] = v[2][0] * (1-weight) + v[3][0] * weight; 
			f2[1] = v[2][1] * (1-weight) + v[3][1] * weight; 
			f2[2] = v[2][2] * (1-weight) + v[3][2] * weight; 

			weight = val[2]/(val[2]-val[6]); 
			f3[0] = v[2][0] * (1-weight) + v[6][0] * weight; 
			f3[1] = v[2][1] * (1-weight) + v[6][1] * weight; 
			f3[2] = v[2][2] * (1-weight) + v[6][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);   geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[4]/(val[4]-val[1]); 

			f4[0] = v[4][0] * (1-weight) + v[1][0] * weight; 
			f4[1] = v[4][1] * (1-weight) + v[1][1] * weight; 
			f4[2] = v[4][2] * (1-weight) + v[1][2] * weight; 

			weight = val[4]/(val[4]-val[3]); 
			f5[0] = v[4][0] * (1-weight) + v[3][0] * weight; 
			f5[1] = v[4][1] * (1-weight) + v[3][1] * weight; 
			f5[2] = v[4][2] * (1-weight) + v[3][2] * weight; 

			weight = val[4]/(val[4]-val[8]); 
			f6[0] = v[4][0] * (1-weight) + v[8][0] * weight; 
			f6[1] = v[4][1] * (1-weight) + v[8][1] * weight; 
			f6[2] = v[4][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f4, f5, f6, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(normal));

			/*	
			if( compute_normal(f4, f5, f6, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f5);   geom->vertex3fv(f6); 
			geom->normal3fv( normal); 
			// 4,5,6
			}
			*/

			weight = val[7]/(val[7]-val[8]); 

			f7[0] = v[7][0] * (1-weight) + v[8][0] * weight; 
			f7[1] = v[7][1] * (1-weight) + v[8][1] * weight; 
			f7[2] = v[7][2] * (1-weight) + v[8][2] * weight; 

			weight = val[7]/(val[7]-val[6]); 
			f8[0] = v[7][0] * (1-weight) + v[6][0] * weight; 
			f8[1] = v[7][1] * (1-weight) + v[6][1] * weight; 
			f8[2] = v[7][2] * (1-weight) + v[6][2] * weight; 

			weight = val[7]/(val[7]-val[3]); 
			f9[0] = v[7][0] * (1-weight) + v[3][0] * weight; 
			f9[1] = v[7][1] * (1-weight) + v[3][1] * weight; 
			f9[2] = v[7][2] * (1-weight) + v[3][2] * weight; 

			compute_normal(f7, f8, f9, normal); 
			vTriangles.push_back(new VECTOR3(f7));
			vTriangles.push_back(new VECTOR3(f8));
			vTriangles.push_back(new VECTOR3(f9));
			vTriangles.push_back(new VECTOR3(normal));

			/*	
			if( compute_normal(f7, f8, f9, normal) != -1 ){
			geom->vertex3fv(f7);
			geom->vertex3fv(f8);   geom->vertex3fv(f9); 
			geom->normal3fv( normal); 
			// 7,8,9
			}
			*/
		}
		no_faces = 3; 
		break;
	case 8:
		{
			double weight = val[1]/(val[1]-val[4]); 
			f1[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			weight = val[2]/(val[2]-val[3]); 
			f2[0] = v[2][0] * (1-weight) + v[3][0] * weight; 
			f2[1] = v[2][1] * (1-weight) + v[3][1] * weight; 
			f2[2] = v[2][2] * (1-weight) + v[3][2] * weight; 

			weight = val[6]/(val[6]-val[7]); 
			f3[0] = v[6][0] * (1-weight) + v[7][0] * weight; 
			f3[1] = v[6][1] * (1-weight) + v[7][1] * weight; 
			f3[2] = v[6][2] * (1-weight) + v[7][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[5]/(val[5]-val[8]); 
			f4[0] = v[5][0] * (1-weight) + v[8][0] * weight; 
			f4[1] = v[5][1] * (1-weight) + v[8][1] * weight; 
			f4[2] = v[5][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f4, f1, f3, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f4, f1, f3, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f1);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 4,1,3
			}
			*/
		}
		no_faces = 2; 
		break;
	case 9:
		{
			double weight = val[1]/(val[1]-val[2]); 
			f1[0] = v[1][0] * (1-weight) + v[2][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[2][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[2][2] * weight; 

			weight = val[6]/(val[6]-val[2]); 
			f2[0] = v[6][0] * (1-weight) + v[2][0] * weight; 
			f2[1] = v[6][1] * (1-weight) + v[2][1] * weight; 
			f2[2] = v[6][2] * (1-weight) + v[2][2] * weight; 

			weight = val[6]/(val[6]-val[7]); 
			f3[0] = v[6][0] * (1-weight) + v[7][0] * weight; 
			f3[1] = v[6][1] * (1-weight) + v[7][1] * weight; 
			f3[2] = v[6][2] * (1-weight) + v[7][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal);
			// 1,2,3
			}
			*/

			weight = val[8]/(val[8]-val[7]); 

			f4[0] = v[8][0] * (1-weight) + v[7][0] * weight; 
			f4[1] = v[8][1] * (1-weight) + v[7][1] * weight; 
			f4[2] = v[8][2] * (1-weight) + v[7][2] * weight; 

			compute_normal(f1, f3, f4, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f3, f4, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f3);
			geom->vertex3fv(f4); 
			geom->normal3fv( normal); 
			// 1,3,4
			}
			*/

			weight = val[1]/(val[1]-val[4]); 
			f5[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f5[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f5[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f1, f4, f5, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f4, f5, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f4);
			geom->vertex3fv(f5); 
			geom->normal3fv( normal); 
			// 1,4,5
			}
			*/

			weight = val[8]/(val[8]-val[4]); 
			f6[0] = v[8][0] * (1-weight) + v[4][0] * weight; 
			f6[1] = v[8][1] * (1-weight) + v[4][1] * weight; 
			f6[2] = v[8][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f5, f4, f6, normal); 
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f5, f4, f6, normal) != -1 ){
			geom->vertex3fv(f5);
			geom->vertex3fv(f4);
			geom->vertex3fv(f6); 
			geom->normal3fv( normal); 
			// 5,4,6
			}
			*/
		}
		no_faces = 4; 
		break;
	case 10:
		{
			double weight = val[1]/(val[1]-val[2]); 
			f1[0] = v[1][0] * (1-weight) + v[2][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[2][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[2][2] * weight; 

			weight = val[4]/(val[4]-val[3]); 
			f2[0] = v[4][0] * (1-weight) + v[3][0] * weight; 
			f2[1] = v[4][1] * (1-weight) + v[3][1] * weight; 
			f2[2] = v[4][2] * (1-weight) + v[3][2] * weight; 

			weight = val[1]/(val[1]-val[5]); 
			f3[0] = v[1][0] * (1-weight) + v[5][0] * weight; 
			f3[1] = v[1][1] * (1-weight) + v[5][1] * weight; 
			f3[2] = v[1][2] * (1-weight) + v[5][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			} 
			*/

			weight = val[4]/(val[4]-val[8]); 
			f4[0] = v[4][0] * (1-weight) + v[8][0] * weight; 
			f4[1] = v[4][1] * (1-weight) + v[8][1] * weight; 
			f4[2] = v[4][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f2, f4, f3, normal); 
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f2, f3, f4, normal) != -1 ){
			geom->vertex3fv(f2);
			geom->vertex3fv(f4);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 2,4,3
			}
			*/

			weight = val[6]/(val[6]-val[2]); 
			f5[0] = v[6][0] * (1-weight) + v[2][0] * weight; 
			f5[1] = v[6][1] * (1-weight) + v[2][1] * weight; 
			f5[2] = v[6][2] * (1-weight) + v[2][2] * weight;  

			weight = val[6]/(val[6]-val[5]); 
			f6[0] = v[6][0] * (1-weight) + v[5][0] * weight; 
			f6[1] = v[6][1] * (1-weight) + v[5][1] * weight; 
			f6[2] = v[6][2] * (1-weight) + v[5][2] * weight;  

			weight = val[7]/(val[7]-val[3]); 
			f7[0] = v[7][0] * (1-weight) + v[3][0] * weight; 
			f7[1] = v[7][1] * (1-weight) + v[3][1] * weight; 
			f7[2] = v[7][2] * (1-weight) + v[3][2] * weight;  

			compute_normal(f5, f6, f7, normal); 
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(f7));
			vTriangles.push_back(new VECTOR3(normal));

			/*
			if( compute_normal(f5, f6, f7, normal) != -1 ){
			geom->vertex3fv(f5);
			geom->vertex3fv(f6);
			geom->vertex3fv(f7); 
			geom->normal3fv( normal); 
			// 5,6,7
			}
			*/

			weight = val[7]/(val[7]-val[8]); 
			f8[0] = v[7][0] * (1-weight) + v[8][0] * weight; 
			f8[1] = v[7][1] * (1-weight) + v[8][1] * weight; 
			f8[2] = v[7][2] * (1-weight) + v[8][2] * weight;  

			compute_normal(f2, f8, f3, normal); 
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f8));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f2, f8, f3, normal) != -1 ){
			geom->vertex3fv(f2);
			geom->vertex3fv(f8);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 2,8,3
			}
			*/
		}
		no_faces = 4; 
		break;
	case 11:
		{
			double weight = val[1]/(val[1]-val[2]); 
			f1[0] = v[1][0] * (1-weight) + v[2][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[2][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[2][2] * weight; 

			weight = val[6]/(val[6]-val[2]); 
			f2[0] = v[6][0] * (1-weight) + v[2][0] * weight; 
			f2[1] = v[6][1] * (1-weight) + v[2][1] * weight; 
			f2[2] = v[6][2] * (1-weight) + v[2][2] * weight; 

			weight = val[7]/(val[7]-val[3]); 
			f3[0] = v[7][0] * (1-weight) + v[3][0] * weight; 
			f3[1] = v[7][1] * (1-weight) + v[3][1] * weight; 
			f3[2] = v[7][2] * (1-weight) + v[3][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[5]/(val[5]-val[8]); 
			f4[0] = v[5][0] * (1-weight) + v[8][0] * weight; 
			f4[1] = v[5][1] * (1-weight) + v[8][1] * weight; 
			f4[2] = v[5][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f1, f3, f4, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f3, f4, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f3);
			geom->vertex3fv(f4); 
			geom->normal3fv( normal); 
			// 1,3,4
			}
			*/

			weight = val[1]/(val[1]-val[4]); 
			f5[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f5[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f5[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f1, f4, f5, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f4, f5, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f4);
			geom->vertex3fv(f5); 
			geom->normal3fv( normal); 
			// 1,4,5
			}
			*/

			weight = val[7]/(val[7]-val[8]); 
			f6[0] = v[7][0] * (1-weight) + v[8][0] * weight; 
			f6[1] = v[7][1] * (1-weight) + v[8][1] * weight; 
			f6[2] = v[7][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f4, f3, f6, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f4, f3, f6, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f3);
			geom->vertex3fv(f6); 
			geom->normal3fv( normal); 
			// 4,3,6
			}
			*/
		}
		no_faces = 4; 
		break;
	case 12:
		{
			double weight = val[2]/(val[2]-val[1]); 
			f1[0] = v[2][0] * (1-weight) + v[1][0] * weight; 
			f1[1] = v[2][1] * (1-weight) + v[1][1] * weight; 
			f1[2] = v[2][2] * (1-weight) + v[1][2] * weight; 

			weight = val[2]/(val[2]-val[3]); 
			f2[0] = v[2][0] * (1-weight) + v[3][0] * weight; 
			f2[1] = v[2][1] * (1-weight) + v[3][1] * weight; 
			f2[2] = v[2][2] * (1-weight) + v[3][2] * weight; 

			weight = val[5]/(val[5]-val[1]); 
			f3[0] = v[5][0] * (1-weight) + v[1][0] * weight; 
			f3[1] = v[5][1] * (1-weight) + v[1][1] * weight; 
			f3[2] = v[5][2] * (1-weight) + v[1][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal);
			// 1,2,3
			}
			*/

			weight = val[5]/(val[5]-val[8]); 
			f4[0] = v[5][0] * (1-weight) + v[8][0] * weight; 
			f4[1] = v[5][1] * (1-weight) + v[8][1] * weight; 
			f4[2] = v[5][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f3, f2, f4, normal); 
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f3, f2, f4, normal) != -1 ){
			geom->vertex3fv(f3);
			geom->vertex3fv(f2);
			geom->vertex3fv(f4); 
			geom->normal3fv( normal); 
			// 3,2,4
			}
			*/

			weight = val[6]/(val[6]-val[7]); 
			f5[0] = v[6][0] * (1-weight) + v[7][0] * weight; 
			f5[1] = v[6][1] * (1-weight) + v[7][1] * weight; 
			f5[2] = v[6][2] * (1-weight) + v[7][2] * weight; 

			compute_normal(f4, f2, f5, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f4, f2, f5, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f2);
			geom->vertex3fv(f5); 
			geom->normal3fv( normal); 
			// 4,2,5
			}
			*/

			weight = val[4]/(val[4]-val[1]); 
			f6[0] = v[4][0] * (1-weight) + v[1][0] * weight; 
			f6[1] = v[4][1] * (1-weight) + v[1][1] * weight; 
			f6[2] = v[4][2] * (1-weight) + v[1][2] * weight; 

			weight = val[4]/(val[4]-val[3]); 
			f7[0] = v[4][0] * (1-weight) + v[3][0] * weight; 
			f7[1] = v[4][1] * (1-weight) + v[3][1] * weight; 
			f7[2] = v[4][2] * (1-weight) + v[3][2] * weight; 

			weight = val[4]/(val[4]-val[8]); 
			f8[0] = v[4][0] * (1-weight) + v[8][0] * weight; 
			f8[1] = v[4][1] * (1-weight) + v[8][1] * weight; 
			f8[2] = v[4][2] * (1-weight) + v[8][2] * weight; 

			compute_normal(f6, f7, f8, normal); 
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(f7));
			vTriangles.push_back(new VECTOR3(f8));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f6, f7, f8, normal) != -1 ){
			geom->vertex3fv(f6);
			geom->vertex3fv(f7);
			geom->vertex3fv(f8); 
			geom->normal3fv( normal); 
			// 6,7,8
			}
			*/
		}
		no_faces = 4; 
		break;
	case 13:
		{
			double weight = val[1]/(val[1]-val[2]); 
			f1[0] = v[1][0] * (1-weight) + v[2][0] * weight; 
			f1[1] = v[1][1] * (1-weight) + v[2][1] * weight; 
			f1[2] = v[1][2] * (1-weight) + v[2][2] * weight; 

			weight = val[1]/(val[1]-val[5]); 
			f2[0] = v[1][0] * (1-weight) + v[5][0] * weight; 
			f2[1] = v[1][1] * (1-weight) + v[5][1] * weight; 
			f2[2] = v[1][2] * (1-weight) + v[5][2] * weight; 

			weight = val[1]/(val[1]-val[4]); 
			f3[0] = v[1][0] * (1-weight) + v[4][0] * weight; 
			f3[1] = v[1][1] * (1-weight) + v[4][1] * weight; 
			f3[2] = v[1][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1);
			geom->vertex3fv(f2);
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[3]/(val[3]-val[2]); 
			f4[0] = v[3][0] * (1-weight) + v[2][0] * weight; 
			f4[1] = v[3][1] * (1-weight) + v[2][1] * weight; 
			f4[2] = v[3][2] * (1-weight) + v[2][2] * weight; 

			weight = val[3]/(val[3]-val[7]); 
			f5[0] = v[3][0] * (1-weight) + v[7][0] * weight; 
			f5[1] = v[3][1] * (1-weight) + v[7][1] * weight; 
			f5[2] = v[3][2] * (1-weight) + v[7][2] * weight; 

			weight = val[3]/(val[3]-val[4]); 
			f6[0] = v[3][0] * (1-weight) + v[4][0] * weight; 
			f6[1] = v[3][1] * (1-weight) + v[4][1] * weight; 
			f6[2] = v[3][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f4, f5, f6, normal); 
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f4, f5, f6, normal) != -1 ){
			geom->vertex3fv(f4);
			geom->vertex3fv(f5);
			geom->vertex3fv(f6); 
			geom->normal3fv( normal); 
			// 4,5,6
			}
			*/

			weight = val[6]/(val[6]-val[2]); 
			f7[0] = v[6][0] * (1-weight) + v[2][0] * weight; 
			f7[1] = v[6][1] * (1-weight) + v[2][1] * weight; 
			f7[2] = v[6][2] * (1-weight) + v[2][2] * weight; 

			weight = val[6]/(val[6]-val[7]); 
			f8[0] = v[6][0] * (1-weight) + v[7][0] * weight; 
			f8[1] = v[6][1] * (1-weight) + v[7][1] * weight; 
			f8[2] = v[6][2] * (1-weight) + v[7][2] * weight; 

			weight = val[6]/(val[6]-val[5]); 
			f9[0] = v[6][0] * (1-weight) + v[5][0] * weight; 
			f9[1] = v[6][1] * (1-weight) + v[5][1] * weight; 
			f9[2] = v[6][2] * (1-weight) + v[5][2] * weight; 

			compute_normal(f7, f8, f9, normal); 
			vTriangles.push_back(new VECTOR3(f7));
			vTriangles.push_back(new VECTOR3(f8));
			vTriangles.push_back(new VECTOR3(f9));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f7, f8, f9, normal) != -1 ){
			geom->vertex3fv(f7);
			geom->vertex3fv(f8);
			geom->vertex3fv(f9); 
			geom->normal3fv( normal); 
			// 7,8,9
			}
			*/

			weight = val[8]/(val[8]-val[5]); 
			f10[0] = v[8][0] * (1-weight) + v[5][0] * weight; 
			f10[1] = v[8][1] * (1-weight) + v[5][1] * weight; 
			f10[2] = v[8][2] * (1-weight) + v[5][2] * weight; 

			weight = val[8]/(val[8]-val[7]); 
			f11[0] = v[8][0] * (1-weight) + v[7][0] * weight; 
			f11[1] = v[8][1] * (1-weight) + v[7][1] * weight; 
			f11[2] = v[8][2] * (1-weight) + v[7][2] * weight; 

			weight = val[8]/(val[8]-val[4]); 
			f12[0] = v[8][0] * (1-weight) + v[4][0] * weight; 
			f12[1] = v[8][1] * (1-weight) + v[4][1] * weight; 
			f12[2] = v[8][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f10, f11, f12, normal); 
			vTriangles.push_back(new VECTOR3(f10));
			vTriangles.push_back(new VECTOR3(f11));
			vTriangles.push_back(new VECTOR3(f12));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f10, f11, f12, normal) != -1 ){
			geom->vertex3fv(f10);
			geom->vertex3fv(f11);
			geom->vertex3fv(f12); 
			geom->normal3fv( normal);
			// 10,11,12
			}
			*/
		}
		no_faces = 4; 
		break;
	case 14:
		{
			double weight = val[2]/(val[2]-val[1]); 
			f1[0] = v[2][0] * (1-weight) + v[1][0] * weight; 
			f1[1] = v[2][1] * (1-weight) + v[1][1] * weight; 
			f1[2] = v[2][2] * (1-weight) + v[1][2] * weight; 

			weight = val[2]/(val[2]-val[3]); 
			f2[0] = v[2][0] * (1-weight) + v[3][0] * weight; 
			f2[1] = v[2][1] * (1-weight) + v[3][1] * weight; 
			f2[2] = v[2][2] * (1-weight) + v[3][2] * weight; 

			weight = val[6]/(val[6]-val[7]); 
			f3[0] = v[6][0] * (1-weight) + v[7][0] * weight; 
			f3[1] = v[6][1] * (1-weight) + v[7][1] * weight; 
			f3[2] = v[6][2] * (1-weight) + v[7][2] * weight; 

			compute_normal(f1, f2, f3, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f2));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f2, f3, normal) != -1 ){
			geom->vertex3fv(f1); 
			geom->vertex3fv(f2); 
			geom->vertex3fv(f3); 
			geom->normal3fv( normal); 
			// 1,2,3
			}
			*/

			weight = val[8]/(val[8]-val[4]); 
			f4[0] = v[8][0] * (1-weight) + v[4][0] * weight; 
			f4[1] = v[8][1] * (1-weight) + v[4][1] * weight; 
			f4[2] = v[8][2] * (1-weight) + v[4][2] * weight; 

			compute_normal(f1, f3, f4, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f3, f4, normal) != -1 ){
			geom->vertex3fv(f1); 
			geom->vertex3fv(f3); 
			geom->vertex3fv(f4); 
			geom->normal3fv( normal); 
			// 1,3,4
			}
			*/

			weight = val[5]/(val[5]-val[1]); 
			f5[0] = v[5][0] * (1-weight) + v[1][0] * weight; 
			f5[1] = v[5][1] * (1-weight) + v[1][1] * weight; 
			f5[2] = v[5][2] * (1-weight) + v[1][2] * weight; 

			compute_normal(f1, f4, f5, normal); 
			vTriangles.push_back(new VECTOR3(f1));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(f5));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f1, f4, f5, normal) != -1 ){
			geom->vertex3fv(f1); 
			geom->vertex3fv(f4); 
			geom->vertex3fv(f5); 
			geom->normal3fv( normal); 
			// 1,4,5
			}
			*/

			weight = val[8]/(val[8]-val[7]); 
			f6[0] = v[8][0] * (1-weight) + v[7][0] * weight; 
			f6[1] = v[8][1] * (1-weight) + v[7][1] * weight; 
			f6[2] = v[8][2] * (1-weight) + v[7][2] * weight; 

			compute_normal(f3, f6, f4, normal); 
			vTriangles.push_back(new VECTOR3(f3));
			vTriangles.push_back(new VECTOR3(f6));
			vTriangles.push_back(new VECTOR3(f4));
			vTriangles.push_back(new VECTOR3(normal));
			/*
			if( compute_normal(f3, f6, f4, normal) != -1 ){
			geom->vertex3fv(f3); 
			geom->vertex3fv(f6); 
			geom->vertex3fv(f4); 
			geom->normal3fv( normal); 
			// 3,6,4
			}
			*/
		}
		no_faces = 4; 
		break;
	default:
		break;
	}

	// if (wcase!=0) glEnd(); 
	//glEnd(); 
	return(no_faces); 
}
