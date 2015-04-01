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

class BoxObject : public Object
{
 private:
  Point3 pmin;
  Point3 pmax;
  int xdim, ydim, zdim, _size;

  vector<Point3> CornerPos;
  vector<Point3> Gradients;
  unsigned short* us_dataPoints;

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
    std::cout << "Calculating corner positions for "<<_size<<" points ..." <<std::endl;

    CornerPos.resize( this->_size );

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
		CornerPos[getIndex(x,y,z)] = Point3(newX, newY, newZ);
	      }
	  }
      }

  }

  void SetData(unsigned short* theData)
  {
    this->us_dataPoints = theData;
    return;
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
       for(int i = 0; i < 36; i+=3){
	 if( Box_IntersectTriangle ( ray, tempHInfo, corners[indices[i]], corners[indices[i+1]], corners[indices[i+2]] )){
	   hitZ.push_back(tempHInfo);
	   boxHit = true;
	   tempHInfo.Init();
	   tempHInfo.node = hInfo.node;
	   tempHInfo.z = hInfo.z;
	 }
       }
         
       if(boxHit) {
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
	 hInfo.shade=RayMarch(ray, hInfo, hitZ[st], hitZ[end], noData);
	 /*hInfo.shade = Color(125,125,0);
	   hInfo.volume = true;
	   hInfo.z = hitZ[st].z;
	   hInfo.N = hitZ[st].N;
	   hInfo.p = hitZ[st].p;
	   return boxHit;*/
	 if(!noData) {
	   //hInfo.volume = true;//false;
	   hInfo.front =  true;
	 }
	 else boxHit = false; //no data (or) iso surface rendering
       }
    return boxHit;
  }
  Color RayMarch(const Ray& ray, HitInfo &hInfo, HitInfo& start, HitInfo& end, bool &noData) const
  {
    Color shade = Color(125,125,0);
    float length = (float)end.z - start.z;
    float dist = sqrt(std::pow((float)end.p.x - start.p.x, 2.0f) + pow((float)end.p.y - start.p.y, 2.0f) + pow((float)end.p.z - start.p.z,  2.0f));

    /* //debug stuff
       shade = Color (dist/10, dist/10, 0.0);
       noData = false;
       hInfo.p = start.p;
       hInfo.z = start.z;
       hInfo.renderIsoSurface = true;
       hInfo.volume = true;
       return shade;
    */

    /** Binning the sample point */
    Ray r;
    r.p = start.p;
    r.dir = end.p - start.p;
    r.Normalize();
    float dt = (float)1.0 / (xdim * 2.0); 
    float t = 0;
    Point3 sample;

    while ( t < dist ) {
      sample = start.p + t * r.dir;
      int cx=-1,cy=-1,cz=-1;

      for(int x=0; x < xdim -1; x++){
	if(CornerPos[getIndex(x,0,0)].x < sample.x && CornerPos[getIndex(x+1,0,0)].x >= sample.x){
	  cx = x;
	  break;
	}
      }
	
      for(int y=0; y < ydim-1; y++){
	if(CornerPos[getIndex(0,y,0)].y < sample.y && CornerPos[getIndex(0,y+1,0)].y >= sample.y){
	  cy = y;
	  break;
	}
      }

      for(int z=0; z < zdim-1; z++){
	if(CornerPos[getIndex(0,0,z)].z < sample.z && CornerPos[getIndex(0,0,z+1)].z >= sample.z){
	  cz = z;
	  break;
	}
      }
	
      if(cx >=0 && cy >=0 && cz >=0 ){
	float data_tot = TrilinearInterpolate(sample,cx,cy,cz);
	float THRESHOLD = 13500;
	if(data_tot > THRESHOLD){
	  noData = false;
	  hInfo.p = sample;
	  hInfo.N = EstimateGradient(cx,cy,cz);
	  hInfo.z += t;
	  hInfo.renderIsoSurface = true;
	  hInfo.volume = false;
	  return shade;
	}
      }	
      t += dt;
	
    }//while loop
    
    noData = true; //didnt find any data that meets the condition
    return shade;
    
  }
  inline int getIndex(int x, int y, int z) const
  {
    return ( z * ydim + y ) * xdim + x;
  }
  float TrilinearInterpolate(Point3 samplePt, int x, int y, int z) const
  {

    float data = 0;
    float xd, yd, zd;
    float c01, c23, c45, c67;
    float c0,c1;
    float val;
    Point3 p1 =  CornerPos[getIndex(x,y,z)];
    Point3 p2 = CornerPos[getIndex(x+1,y+1,z+1)];
    float* v = new float[8];

    v[0] = us_dataPoints[getIndex(x,y,z)];
    v[1] = us_dataPoints[getIndex(x+1,y,z)];
    v[2] = us_dataPoints[getIndex(x,y+1,z)];
    v[3] = us_dataPoints[getIndex(x+1,y+1,z)];
    v[4] = us_dataPoints[getIndex(x+1,y,z+1)];
    v[5] = us_dataPoints[getIndex(x,y,z+1)];
    v[6] = us_dataPoints[getIndex(x,y+1,z+1)];
    v[7] = us_dataPoints[getIndex(x+1,y+1,z+1)];

    xd = (samplePt.x - p1.x)/(p2.x - p1.x);
    yd = (samplePt.y - p1.y)/(p2.y - p1.y);
    zd = (samplePt.z - p1.z)/(p2.z - p1.z);
    c01 = v[0] * (1 - xd) + v[1] * xd;
    c23 = v[2] * (1 - xd) + v[3] * xd;
    c45 = v[4] * (1 - xd) + v[5] * xd;
    c67 = v[6] * (1 - xd) + v[7] * xd;
    c0 = c01 * (1 - yd) + c45 * yd;
    c1 = c23 * (1 - yd) + c67 * yd;
    val = c0 * (1 - zd) + c1 * zd;

    delete[] v;
    return val;
  }
  void CalculateGradients()
  {
    std::cout<<"Calculating Gradient differences ..."<<std::endl;

    Gradients.resize(this->_size);
    for (int z = 0; z < zdim; z++)
      {
	for(int y=0;y<ydim;y++)
	  {
	    for(int x=0; x<xdim;x++)
	      {
		float xn, xp, yn, yp, zn, zp;
		if(x > 0) xp = us_dataPoints[getIndex(x-1,y,z)];//DataPoints[x-1][y][z];
		else xp=0;
		if(x < xdim-1) xn = us_dataPoints[getIndex(x+1,y,z)];//DataPoints[x+1][y][z];
		else xn = 0;
		if(y > 0) yp = us_dataPoints[getIndex(x,y-1,z)];//DataPoints[x][y-1][z];
		else  yp = 0;
		if(y < ydim-1) yn = us_dataPoints[getIndex(x,y+1,z)];//DataPoints[x][y+1][z];
		else yn = 0;
		if(z > 0) zp = us_dataPoints[getIndex(x,y,z-1)];//DataPoints[x][y][z-1];
		else zp=0;
		if(z < zdim-1) zn = us_dataPoints[getIndex(x,y,z+1)];//DataPoints[x][y][z+1];
		else zn = 0;
		Point3 N =  Point3(xn - xp, yn - yp, zn - zp);
		Gradients[getIndex(x,y,z)] = N;
	      }
	  }
      }

  }

  Point3 EstimateGradient(int x, int y, int z) const
  {
    /* maybe i can calculate this avg while calculating gradients; maybe needto trilinear interpolate gradient
       as well */
    Point3 grad = Gradients[getIndex(x,y,z)];
    grad += Gradients[getIndex(x+1,y,z)];
    grad += Gradients[getIndex(x,y+1,z)];
    grad += Gradients[getIndex(x+1,y+1,z)];
    grad += Gradients[getIndex(x+1,y,z+1)];
    grad += Gradients[getIndex(x,y,z+1)];
    grad += Gradients[getIndex(x,y+1,z+1)];
    grad += Gradients[getIndex(x+1,y+1,z+1)];

    grad.x /= -8.0; grad.y /= -8.0; grad.z /= -8.0;
    if(grad.x != 0 && grad.y!= 0 && grad.z!=0)
      grad.Normalize();
    return grad;
    
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
