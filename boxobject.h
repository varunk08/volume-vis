//Box object for ray-tracing
#ifndef BOXOBJECT_H
#define BOXOBJECT_H

#include <vector>
#include <iostream>
#include <string.h>
#include <cmath>
#include <GLUT/GLUT.h>
#include "objects.h"
#include "scene.h"

using namespace std;
typedef struct
{
  vector <vector <int> > Indices;
} Cell;

class BoxObject : public Object
{
 private:
  Point3 pmin;
  Point3 pmax;
  int xdim, ydim, zdim, _size;
  //Corner points data
  vector< vector<vector <Point3> > > CornerLocations;
  //Data points
  vector< vector<vector <int> > > DataPoints;
  //Cells
  vector < vector <vector <Cell> > > Cells;
 public:
  BoxObject()
    {
      pmin = Point3(-1, -1, -1);
      pmax = Point3(1, 1, 1);
    }
  Box GetBoundBox() const
  {
    return Box(-1,-1,-1,1,1,1);
  }
  void SetDimensions(int xdim, int ydim, int zdim)
  {
    this->xdim = xdim;
    this->ydim = ydim;
    this->zdim = zdim;
    this->_size = xdim * ydim * zdim;
    std::cout << "Setting dimensions-> Size: "<<_size<<std::endl;
    //create cells
    std::cout<<"Filling corner locations ..."<<std::endl;
    CornerLocations.resize(xdim);
    for(int i = 0; i < xdim; i++)
      {
	CornerLocations[i].resize(ydim);
	for(int j = 0; j < ydim; j++)
	  CornerLocations[i][j].resize(zdim);
      }
    float startX = -1.0f;
    float startY = -1.0f;
    float startZ = -1.0f;
    float newX = 0.0f, newY = 0.0f, newZ = 0.0f;
    for(int z = 0; z < zdim; z++)
      {
	newZ = startZ + z * (2.0 / (zdim - 1));
	for(int y = 0; y < ydim; y++)
	  {
	    newY = startY + y * (2.0 / (ydim - 1));
	    for(int x = 0; x < xdim; x++)
	      {
		newX = startX + x * (2.0 / (xdim - 1));
		CornerLocations[x][y][z] = Point3(newX, newY, newZ);
		//		std::cout<<"CornerLocations["<<x<<"]["<<y<<"]["<<z<<"]: "<<CornerLocations[x][y][z].x<<","<<CornerLocations[x][y][z].y<<","<<CornerLocations[x][y][z].z<<std::endl;
	      }
	  }
      }

    Point3 c = CornerLocations[xdim-1][ydim-1][zdim-1];
    std::cout<<"Corner: " <<c.x << " "<<c.y<<" "<<c.z<<" "<<std::endl;


  }
  void SetData(unsigned short* theData)
  {
    std::cout<<"Filling data ..."<<std::endl;
    unsigned int size = xdim * ydim * zdim;
    std::cout<<"Size: "<<size<<std::endl;
    DataPoints.resize(xdim);
    for(int i = 0; i < xdim; i++)
      {
	DataPoints[i].resize(ydim);
	for(int j = 0; j < ydim; j++)
	  DataPoints[i][j].resize(zdim);
      }
    int index = 0;
    for(int z=0; z < zdim; z++)
      {
	for(int y = 0; y < ydim; y++)
	  {
	    for(int x = 0; x < xdim; x++)
	      {
		DataPoints[x][y][z] =(int) theData[index];
		//		std::cout<<"DataPoints["<<x<<"]["<<y<<"]["<<z<<"]: "<<theData[index]<<std::endl;
		index++;
	      }
	  }
	
      }

    //    std::cout<<DataPoints[0][0][0]<<std::endl;

  }
  void CreateCells()
  {
    std::cout<<"Creating cells ..."<<std::endl;
    int ix = 0, iy = 0, iz= 0;
    int n_cells_x = xdim - 1;
    int n_cells_y = ydim - 1;
    int n_cells_z = zdim - 1;
    Cells.resize(n_cells_x);
    for(int i = 0; i < n_cells_x; i++)
      {
	Cells[i].resize(n_cells_y);
	for(int j=0; j<n_cells_y; j++)
	  Cells[i][j].resize(n_cells_z);
      }

    unsigned int ncells = 0; unsigned int totalCells = n_cells_x * n_cells_y * n_cells_z;
    for(int z = 0; z < n_cells_z; z++)
      {
	iy=0;
	for(int y = 0; y < n_cells_y; y++)
	  {
	    ix = 0;
	    for(int x = 0; x < n_cells_x; x++)
	      {
		Cell c;
		c.Indices.resize(8);
		for(int i = 0; i< 8;i++)
		  {
		    c.Indices[i].resize(3);
		  }

	       	c.Indices[0] = {ix,iy,iz};
		c.Indices[1] = {ix+1,iy,iz};
		c.Indices[2] = {ix,iy+1,iz};
		c.Indices[3] = {ix+1,iy+1,iz};
		c.Indices[4] = {ix+1,iy,iz+1};
		c.Indices[5] = {ix,iy,iz+1};
		c.Indices[6] = {ix,iy+1,iz+1};
		c.Indices[7] = {ix+1,iy+1,iz+1};
		Cells[x][y][z] = c;
		ix++;
		//ncells++;
		
		//		std::cout<<"Progress: " <<ncells<<"/"<<totalCells<<"\r";
		//	        std::cout.flush();
	      }
	    iy++;
	  }
	iz++;
      }
    std::cout<<std::endl;
    std::cout<<Cells[0][0].size()<<std::endl;;
    std::cout<<"Test:\n"<<std::endl;
		Cell c = Cells[0][0][0];
				  Point3 corners[8];
				  corners[0] = Point3(CornerLocations[c.Indices[0][0]][c.Indices[0][1]][c.Indices[0][2]]);
				  corners[1] = Point3(CornerLocations[c.Indices[1][0]][c.Indices[1][1]][c.Indices[1][2]]);
				  corners[2] = Point3(CornerLocations[c.Indices[2][0]][c.Indices[2][1]][c.Indices[2][2]]);
				  corners[3] = Point3(CornerLocations[c.Indices[3][0]][c.Indices[3][1]][c.Indices[3][2]]);
				  corners[4] = Point3(CornerLocations[c.Indices[4][0]][c.Indices[4][1]][c.Indices[4][2]]);
				  corners[5] = Point3(CornerLocations[c.Indices[5][0]][c.Indices[5][1]][c.Indices[5][2]]);
				  corners[6] = Point3(CornerLocations[c.Indices[6][0]][c.Indices[6][1]][c.Indices[6][2]]);
				  corners[7] = Point3(CornerLocations[c.Indices[7][0]][c.Indices[7][1]][c.Indices[7][2]]);
				  //       for(int i = 0; i < 8; i++)
				  //        {
				  //	  std::cout<<corners[i].x<<" "<<corners[i].y<<" "<<corners[i].z<<std::endl;
				  //	}
       //       std::cout<<CornerLocations.size()<<std::endl;
       //       std::cout<<CornerLocations[0].size()<<std::endl;
       //       std::cout<<CornerLocations[0][0].size()<<std::endl;
  }
  bool IntersectRay(const Ray &ray, HitInfo &hInfo, int hitSide=HIT_FRONT ) const
  {
    bool boxHit = false;
    //the ray is already in object space.
    Point3 corners[8];
    corners[0] = Point3(-1.0f, -1.0f, -1.0f);
    corners[1] = Point3(1.0f, -1.0f, -1.0f);
    corners[2] = Point3(1.0f, -1.0f, 1.0f);
    corners[3] = Point3(-1, -1, 1);
    corners[4] = Point3(-1, 1, -1);
    corners[5] = Point3(1, 1, -1);
    corners[6] = Point3(1, 1, 1);
    corners[7] = Point3(-1, 1, 1);
    int indices[36] = {
      0,1,2,
      0,2,3,//bottom
      3,0,4,
      3,4,7,//left
      4,5,6,
      4,6,7,//top
      5,2,6,
      5,1,2,//right
      0,5,4,
      0,1,5,//front
      2,3,7,
      2,7,6//back
    };
    HitInfo tempHInfo;
    tempHInfo.Init();
    tempHInfo.node = hInfo.node;
    tempHInfo.z = hInfo.z;
    std::vector<HitInfo> hitZ;
    hitZ.clear();
       for(int i = 0; i < 36; i+=3)
        {
	 if( Box_IntersectTriangle ( ray, tempHInfo, corners[indices[i]], corners[indices[i+1]], corners[indices[i+2]] ))
	   {
	    hitZ.push_back(tempHInfo);
	      boxHit = true;
	      tempHInfo.Init();
	      tempHInfo.node = hInfo.node;
	      tempHInfo.z = hInfo.z;
	    }
	  }
         
       ///std::cout<<"box hit z: "<<hitZ[0].z<<" " <<hitZ[1].z<<std::endl;
       if(boxHit)
	 {
	   //std::cout<<"Box hit "<<hitZ.size()<<std::endl;
	       int st = 0, end = 0;
	       if( hitZ[0].z < hInfo.z && hitZ[1].z < hInfo.z ) {
		 if( hitZ[0].z < hitZ[1].z ) { 
		   st = 0; end = 1;
		 }
		 else {
		   st = 1; end = 0;
		 }
	       }
	       bool noData = false;
	       hInfo.shade=RayMarch(ray, hitZ[st], hitZ[end], noData);
	       if(!noData) {
	       //hInfo.p = hitZ[ci].p;
	       //hInfo.N = hitZ[ci].N;
	       hInfo.volume = true;
	       //hInfo.node = hitZ[ci].node;
	       //hInfo.z = hitZ[ci].z;
	       hInfo.front =  true;
	       }
	       else boxHit = false;
	 }
    //create 12 triangles with the corners of the box
    //do intersection testing with all the triangles and get correct normal for correct shading

    return boxHit;
  }
  Color RayMarch(const Ray& ray, HitInfo& start, HitInfo& end, bool &noData) const
  {
    Color shade = Color(125,125,0);
    float length = (float)end.z - start.z;
    //    std::cout<<start.p.x<<" " <<start.p.y<<" "<<start.p.z<<std::endl;
    //    std::cout<<end.p.x<<" " <<end.p.y<<" "<<end.p.z<<std::endl;
    //    std::cout<<end.z<<" "<<start.z<<" " << std::endl;
    float dist = sqrt(pow(end.p.x - start.p.x, 2.0) + pow(end.p.y - start.p.y, 2.0) + pow(end.p.z - start.p.z,  2.0));
    shade = Color (dist/10, dist/10, 0.0);

    //we know the ray. start the binning process- find which cell the current sample point is in, get the data points and average the color
    Ray r;
    r.p = start.p;
    r.dir = end.p - start.p;
    r.Normalize();
    float dt = dist / (2.0); //testing
    //get the first sample only for now -- testing
    Point3 sample = start.p + dt * r.dir;
    //bin x
    int cx,cy,cz;
	    for(int x =0; x < xdim  -1; x++)
	      {
		Cell c = Cells[x][0][0];
		if( CornerLocations[c.Indices[0][0]][c.Indices[0][1]][c.Indices[0][2]].x < sample.x && CornerLocations[c.Indices[7][0]][c.Indices[7][1]][c.Indices[7][2]].x > sample.x)
		  {
		    cx = x;
		  } 
	      }
	    for(int y =0; y < ydim -1; y++)
	      {
		Cell c = Cells[0][y][0];
		if( CornerLocations[c.Indices[0][0]][c.Indices[0][1]][c.Indices[0][2]].y < sample.y && CornerLocations[c.Indices[7][0]][c.Indices[7][1]][c.Indices[7][2]].y > sample.y)
		  {
		    cy = y;
		  } 

	      }
	    for(int z = 0; z < zdim -1; z++)
	      {
		Cell c = Cells[0][0][z];
		if( CornerLocations[c.Indices[0][0]][c.Indices[0][1]][c.Indices[0][2]].z < sample.z && CornerLocations[c.Indices[7][0]][c.Indices[7][1]][c.Indices[7][2]].z > sample.z)
		  {
		    cz = z;
		  } 
	      }
    Cell tc = Cells[cx][cy][cz];
    //average the datapoints
    float data_tot = 0;
    //for(int i =0; i<8; i++)
      {
    std::cout<<cx<<","<<cy<<","<<cz<< " Data: "<<std::endl;
	data_tot = DataPoints[tc.Indices[0][0]][tc.Indices[0][1]][tc.Indices[0][2]];
      }

      //data_tot /= 8;

    data_tot /= 1000;
    //    if( data_tot <= 500 ) noData = true;
    noData = false;
    shade = Color(data_tot, data_tot, data_tot);
    return shade;

  } 

  bool Box_IntersectTriangle(const Ray &ray, HitInfo &hInfo, Point3 p1, Point3 p2, Point3 p3, bool raymarch=false) const
  {
    //Based on Shirley's book
    Point3 A, B, C;
    A = p1; B = p2; C = p3;
    float a = A.x - B.x;
    float b = A.y - B.y;
    float c = A.z - B.z;

    float d = A.x - C.x;
    float e = A.y - C.y;
    float f = A.z - C.z;

    float g = ray.dir.x;
    float h = ray.dir.y;
    float i = ray.dir.z;

    float j = A.x - ray.p.x;
    float k = A.y - ray.p.y;
    float l = A.z - ray.p.z;

    float eimhf = e * i - h * f;
    float gfmdi = g * f - d * i;
    float dhmeg = d * h - e * g;
    float akmjb = a * k - j * b;
    float jcmal = j * c - a * l;
    float blmkc = b * l - k * c;

    float M = a * (eimhf) + b * (gfmdi) + c * (dhmeg);
    float t = -(f * akmjb + e * jcmal + d * blmkc)/M;

      if( t > BIGFLOAT || t > hInfo.z)
	return false;

    float gamma = (i*(akmjb) + h * ( jcmal ) + g * ( blmkc ))/M;

    if( gamma  < 0 || gamma > 1)
      return false;

    float beta = (j * (eimhf) + k * gfmdi + l * dhmeg )/M;

    if(beta < 0 || beta > 1 - gamma)
      return false;

    hInfo.z = t;
    hInfo.p = ray.p + t * ray.dir;
    Point3 N = (C - A).Cross(B - A);
    N.Normalize();
    hInfo.N = N;
    hInfo.volume = true;

    
    return true;
    /*
    Point3 A, B, C, N, P, D, H;
    float bias = 1e-6f;
    float t, a, a1, a2, alpha, beta;
    Point2 A2d, B2d, C2d, H2d;


    P = ray.p;
    D = ray.dir;
    Point3 nTemp = (B-A).Cross(C-A);
    nTemp.Normalize();
    N = nTemp;
    
    if(N.Dot(P-A) < bias) return false;
    if(N.Dot(D) == 0) return false;
    t = N.Dot(P-A)/ -(N.Dot(D));

  if( t < hInfo.z && t >= bias && t < BIGFLOAT ) 
    {
      H = P + t * D;
      int projAxis = 0;
            if(abs(H.x) > abs(H.y)){
                if(abs(H.x) > abs(H.z)){
                    projAxis = 0;
                }
                else{
                    projAxis = 2;
                }
            }
            else if(abs(H.y) > abs(H.z)){
                projAxis = 1;
            }
            else{
                projAxis = 2;
            }
            //cout<<projAxis<<endl;
            switch (projAxis) {
                case 0: 
                    A2d = Point2(A.y, A.z);
                    B2d = Point2(B.y, B.z);
                    C2d = Point2(C.y, C.z);
                    H2d = Point2(H.y, H.z);
                    break;
                case 1:// project onto the y-axis
                    A2d = Point2(A.x, A.z);
                    B2d = Point2(B.x, B.z);
                    C2d = Point2(C.x, C.z);
                    H2d = Point2(H.x, H.z);
                    break;
                case 2:// project onto the z-axis
                    A2d = Point2(A.x, A.y);
                    B2d = Point2(B.x, B.y);
                    C2d = Point2(C.x, C.y);
                    H2d = Point2(H.x, H.y);
                default:
                    break;
            }
            
            a = (B2d - A2d).Cross(C2d - A2d)/ 2.0;
            a1 = (B2d - A2d).Cross(H2d - A2d)/ 2.0;
            a2 = (H2d - A2d).Cross(C2d - A2d)/ 2.0;
            alpha = a1/a;
            beta = a2/a;
            if(alpha < -bias || beta < -bias || alpha + beta > 1+bias){
                return false; // point H is outside triangle
            }
            float b0 = 1 - beta - alpha;
            float b1 = beta;
            float b2 = alpha;
            H = b0 * A + b1 * B + b2 * C;

	    hInfo.front = true;
	    hInfo.p = H;
	    hInfo.N = N;
	    hInfo.z = t;

	    return true;
	    
    }
  return false;
    */
  }
bool IsInside(const Point3 &p) const { for ( int i=0; i<3; i++ ) if ( pmin[i] > p[i] || pmax[i] < p[i] ) return false; return true; }
  void ViewportDisplay() const
  {
    /*        static GLUquadric *q = NULL;
        if ( q == NULL ) q = gluNewQuadric();
        gluSphere(q,1,50,50);

    */
    //Multi-colored side - FRONT
    glBegin(GL_POLYGON);
 
    glColor3f( 1.0, 0.0, 0.0 );     glVertex3f(  1.0, -1.0, -1.0 );      // P1 is red
    glColor3f( 0.0, 1.0, 0.0 );     glVertex3f(  1.0,  1.0, -1.0 );      // P2 is green
    glColor3f( 0.0, 0.0, 1.0 );     glVertex3f( -1.0,  1.0, -1.0 );      // P3 is blue
    glColor3f( 1.0, 0.0, 1.0 );     glVertex3f( -1.0, -1.0, -1.0 );      // P4 is purple
 
    glEnd();
 
    // White side - BACK
    glBegin(GL_POLYGON);
    glColor3f(   1.0,  1.0, 1.0 );
    glVertex3f(  1.0, -1.0, 1.0 );
    glVertex3f(  1.0,  1.0, 1.0 );
    glVertex3f( -1.0,  1.0, 1.0 );
    glVertex3f( -1.0, -1.0, 1.0 );
    glEnd();
 
    // Purple side - RIGHT
    glBegin(GL_POLYGON);
    glColor3f(  1.0,  0.0,  1.0 );
    glVertex3f( 1.0, -1.0, -1.0 );
    glVertex3f( 1.0,  1.0, -1.0 );
    glVertex3f( 1.0,  1.0,  1.0 );
    glVertex3f( 1.0, -1.0,  1.0 );
    glEnd();
 
    // Green side - LEFT
    glBegin(GL_POLYGON);
    glColor3f(   0.0,  1.0,  0.0 );
    glVertex3f( -1.0, -1.0,  1.0 );
    glVertex3f( -1.0,  1.0,  1.0 );
    glVertex3f( -1.0,  1.0, -1.0 );
    glVertex3f( -1.0, -1.0, -1.0 );
    glEnd();
 
    // Blue side - TOP
    glBegin(GL_POLYGON);
    glColor3f(   0.0,  0.0,  1.0 );
    glVertex3f(  1.0,  1.0,  1.0 );
    glVertex3f(  1.0,  1.0, -1.0 );
    glVertex3f( -1.0,  1.0, -1.0 );
    glVertex3f( -1.0,  1.0,  1.0 );
    glEnd();
 
    // Red side - BOTTOM
    glBegin(GL_POLYGON);
    glColor3f(   1.0,  0.0,  0.0 );
    glVertex3f(  1.0, -1.0, -1.0 );
    glVertex3f(  1.0, -1.0,  1.0 );
    glVertex3f( -1.0, -1.0,  1.0 );
    glVertex3f( -1.0, -1.0, -1.0 );
    glEnd();
 
    glFlush();

  }
};
extern BoxObject theBoxObject;
#endif
