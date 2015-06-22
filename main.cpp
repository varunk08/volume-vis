#include "tinyxml/tinyxml.h"
//#include <GLUT/GLUT.h>
//#include "viewport.h"
#include "utils.h"
#include "scene.h"
#include "xmlload.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include "cyColor.h"

#include <unistd.h>
#include <vector>
#include "threadpool.h"
#include "volumedata.h"
#include "boxobject.h"

using namespace std;

#define BIAS_SHADOW 1e-4f
#define BIAS_SHADING 0.001f
#define _USE_MATH_DEFINES
#define NUM_THREADS  8

Camera camera;
Node rootNode;
RenderImage renderImage;

BoxObject theBoxObject;
Sphere theSphere;
Plane thePlane;
MaterialList materials;
LightList lights;
ObjFileList objList;




Color24 white(cyColor(255.0));
Color24 black(cyColor(0.0));

bool TraceNode(const Ray &r,HitInfo &hitInfo,const Node &curNode);
bool RayTrace(HitInfo &hitInfo, Node* node,Ray ray,int PixIndex);
bool RayTrace_2(const Ray &ray, HitInfo &hitInfo);


//! Un-threaded begin render
void BeginRender()
{
  std::cout<<"Rendering ...\n";
    
  float alpha = camera.fov;
  float l = 1.0;
  float h = l * tan(alpha/2.0 *(M_PI/180));
  float aspectRatio = (float)camera.imgWidth/camera.imgHeight;
  float s = aspectRatio * abs(h);
  float dx = (2 * abs(s))/camera.imgWidth;
  float dy = -(2 * abs(h))/camera.imgHeight;
  float dxx = dx/2,dyy=dy/2;
  Point3 K(-s,h,-l);
    
  K.x += dxx;
  K.y += dyy;
    
  for(int i = 0; i< camera.imgHeight; i++){
    for(int j = 0; j<camera.imgWidth; j++){
      std::cout<<"("<<j<<" / "<<camera.imgWidth<<","<<i<<" / "<<camera.imgHeight<<")\r";
      std::cout.flush();
      int PixIndex = i * camera.imgWidth + j;
      bool pixelHit=false;
      K.x += dx;
      Matrix3 RotMat;
      cyPoint3f f = camera.dir;
      f.Normalize();
      cyPoint3f s = f.Cross(camera.up);
      s.Normalize();
      cyPoint3f u = s.Cross(f);
      const float pts[9]={s.x,u.x,-f.x,s.y,u.y,-f.y,s.z,u.z,-f.z};
      RotMat.Set(pts);
    
      Ray r(camera.pos, K);
      r.dir=r.dir*RotMat;
      r.dir.Normalize();
    
      HitInfo hitInfo;
      hitInfo.Init();
           
    
      Color shade(255,255,255);
      if(rootNode.GetNumChild()>0){
	if(RayTrace_2(r, hitInfo)) {
	  pixelHit=true;
	  if(hitInfo.volume) {
	    //shade = hitInfo.node->GetMaterial()->VolumeShade(r, hitInfo, lights, hitInfo.shade);
	    shade = hitInfo.shade;
	    }
	  else{ //might render ISO surface
	    //if(hitInfo.renderIsoSurface) std::cout<<"shading iso"<<std::endl;
	    shade = hitInfo.node->GetMaterial()->Shade(r, hitInfo, lights, 5);
	  }
	}
      
	renderImage.PutPixel(PixIndex, shade, hitInfo.z);
      }
      if(!pixelHit){
	renderImage.PutPixel(PixIndex, black, BIGFLOAT);
      }
    }
    K.x = -s;
    K.x += dxx;
    K.y += dy;
  }
  cout<<"Render Complete.\n"<<endl;
  renderImage.ComputeZBufferImage();
  renderImage.SaveZImage("images/zbuffer.ppm");
  renderImage.SaveImage("images/renderimage.ppm");
    
}

float GenLight::Shadow(Ray ray, float t_max)
{
   // cout<<"Calculating shadow"<<endl;
   //DISABLING SHADOW CALCULATIONS -- VOLUME RENDERING
  /*    float eps = BIAS_SHADOW;
    ray.p = Point3(ray.p.x+eps*ray.dir.x, ray.p.y+eps*ray.dir.y,ray.p.z+eps*ray.dir.z);
    HitInfo hitInfo;
    
    if(RayTrace_2(ray, hitInfo))
    {
        
        if(hitInfo.z>0 && hitInfo.z < t_max){
            
            return 0.0; //occluded
        }
        
	}*/
    return 1.0; //direct
}

bool Box::IntersectRay(const Ray &ray, float t_max) const
{
    Ray r = ray;
    float tmin  =  -t_max;
    float tmax = t_max;
    //if ray is inside box - return true
    if (IsInside(r.p)) return true;
    //get pairs of planes - x , y, z
    // 0:(x_min,y_min,z_min), 1:(x_max,y_min,z_min)
    // 2:(x_min,y_max,z_min), 3:(x_max,y_max,z_min)
    // 4:(x_min,y_min,z_max), 5:(x_max,y_min,z_max)
    // 6:(x_min,y_max,z_max), 7:(x_max,y_max,z_max)
    float xl = Corner(0).x;
    float xh = Corner(3).x;
    float yl = Corner(0).y;
    float yh = Corner(2).y;
    float zl = Corner(0).z;
    float zh = Corner(5).z;
    
    //Check intersection for X planes
    if(r.p.x != 0.0){
        float tx1 = (xl - ray.p.x)/ ray.dir.x;
        float tx2 = (xh - ray.p.x)/ ray.dir.x;
        
        tmin = utils::max((float)tmin,(float) utils::min((float)tx1,tx2));
        tmax = utils::min((float)tmax,(float) utils::max((float)tx1, tx2));
    }
    //Check intersection for Y planes
    if(r.p.y != 0.0){
        float tx1 = (yl - ray.p.y)/ ray.dir.y;
        float tx2 = (yh - ray.p.y)/ ray.dir.y;
        
        tmin = max(tmin, min(tx1,tx2));
        tmax = min(tmax, max(tx1, tx2));
    }
    
    //Check intersection for Z planes
    if(r.p.z != 0.0){
        float tx1 = (zl - ray.p.z)/ ray.dir.z;
        float tx2 = (zh - ray.p.z)/ ray.dir.z;
        
        tmin = max(tmin, min(tx1,tx2));
        tmax = min(tmax, max(tx1, tx2));
    }
    
    return tmax>=tmin;
    
}
Color MtlBlinn::VolumeShade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, Color preCol) const
{
      Color shade;

    Color ambComponent = Color(0,0,0);
    Color diffuse = preCol;
    Color specular = Color(0.7f,0.7f,0.7f);
    Color ambInt = preCol;
    Color allOther = Color(0,0,0);
    float glossiness = 20.0f;    
    Point3 P = Point3(hInfo.p);
    for ( unsigned int i=0; i<lights.size(); i++ ) {
        if(lights[i]->IsAmbient()){
            Color intensity = lights[i]->Illuminate(hInfo.p);
            ambComponent += (ambInt * intensity);
            continue;
        }
        else{
            Point3 L = -lights[i]->Direction(P);
            L.Normalize();
            
            Point3 V = ray.p - P;
            V.Normalize();
            
            Point3 LplusV = L + V;
            Point3 H = (L+V)/LplusV.Length();
            H.Normalize();
            
            float alpha = glossiness;
            Point3 N = hInfo.N;
            float S = H.Dot(N);
            S = pow((float)S,alpha);
            float costheta = L.Dot(N)/(L.Length() * N.Length());
            Color intensity = lights[i]->Illuminate(P);
            allOther += intensity * (costheta>0?costheta:0)  * (diffuse + S * (specular)) ;
        }
        /* finally add inta*cola + intall*costheta*(cold + s* colS)*/
        shade = ambComponent  + allOther;
    }

    return shade;

}
Color MtlBlinn::Shade(const Ray &ray, const HitInfo &hInfo, const LightList &lights, int bounceCount) const{
    float bias = BIAS_SHADING;
    Color shade;
    Color rShade = Color(0,0,0);
    Color tShade = Color(0,0,0);
    const Material *mat;
    mat = hInfo.node->GetMaterial();
    const MtlBlinn* mb =static_cast<const MtlBlinn*>(mat);
//    cout<<"HInfo front: "<<hInfo.front<<endl;
    /* local copy */
    Point3 P;
    P.Set(hInfo.p.x,hInfo.p.y,hInfo.p.z);
    Ray iRay = ray;
    
    Color ambInt = mb->diffuse;
    Color allOther = Color(0,0,0);
    Color diffuse = mb->diffuse;;
    Color ambComponent = Color(0,0,0);
    
    for ( unsigned int i=0; i<lights.size(); i++ ) {
        if(lights[i]->IsAmbient()){
//            cout<<"ambient "<<endl;
            Color intensity = lights[i]->Illuminate(hInfo.p);
            ambComponent += (ambInt * intensity);
            continue;
        }
        else{
//            cout<<"other lighting  "<<endl;
            Point3 L = -lights[i]->Direction(P);
            L.Normalize();
            
            Point3 V = ray.p - P;
            V.Normalize();
            
            Point3 LplusV = L + V;
            Point3 H = (L+V)/LplusV.Length();
            H.Normalize();
            
            float alpha = mb->glossiness;
            Point3 N = hInfo.N;
            float S = H.Dot(N);
            S = pow((float)S,alpha);
            float costheta = L.Dot(N)/(L.Length() * N.Length());
            Color intensity = lights[i]->Illuminate(P);
//            cout<<"costheta "<<endl;
            allOther += intensity * (costheta>0?costheta:0) * (diffuse + S * (mb->specular)) ;
        }
        /* finally add inta*cola + intall*costheta*(cold + s* colS)*/
        shade = ambComponent  + allOther;
    }

    return shade; //testing -- no need reflection or refraction now    
    /* Calculate refraction */
    if(refraction.Grey()>0 && bounceCount>0){
        Color reflShade = Color(0,0,0);
        float R0, Refl = 0.0f, Trans = 0.0f;
        HitInfo temp;
        temp.Init();
        
        Point3 N = hInfo.N;
//        Point3 V = Point3(iRay.p.x -  hInfo.p.x, iRay.p.y - hInfo.p.y, iRay.p.z - hInfo.p.z);
        Point3 V = Point3(hInfo.p.x - iRay.p.x, hInfo.p.y - iRay.p.y, hInfo.p.z - iRay.p.z);
        V.Normalize();
        float n1 = 1, n2 = 1;
        if(hInfo.front){ /* Hitting from outside */
//            temp.front = false;
            n2 = ior;
//            cout<<"outside "<<endl;
        }
        else if(!hInfo.front){ /* Transmission from the inside */
//            temp.front = true;
            n1 = ior;
//            cout<<"intside... "<<endl;
            N = -hInfo.N;
        }
        float ratio_n = n1 / n2;
        
        float costheta_v = -V.Dot(N);        /* refer: http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf */

        float sin2theta_t = ratio_n * ratio_n * (1 - costheta_v * costheta_v);
        Point3 T =   ratio_n * V + (ratio_n * costheta_v - sqrtf(1 - sin2theta_t)) * N ;
//        cout<<ratio_n<<" "<<"cos_v "<<costheta_v<<" sin2theta_t "<<sin2theta_t<<endl;
        Ray tRay = Ray(hInfo.p,T);
        
        //tRay.dir.Normalize();
        tRay.p.x = tRay.p.x + bias *tRay.dir.x; /* add bias */
        tRay.p.y = tRay.p.y + bias *tRay.dir.y;
        tRay.p.z = tRay.p.z + bias *tRay.dir.z;
//        cout<<"B temp front: "<< temp.front<<endl;
        if(sin2theta_t <= 1){
            if(RayTrace_2(tRay, temp)){
//                bounceCount--;
//                cout<<"A temp front: "<< temp.front<<endl;
                tShade  =  temp.node->GetMaterial()->Shade(tRay,temp,lights,bounceCount);
                tShade.r *= exp(-absorption.r * temp.z);
                tShade.g *= exp(-absorption.g * temp.z);
                tShade.b *= exp(-absorption.b * temp.z);
//                shade = tShade; /* remove later */
//                return shade;
               
                
                /* Calculate Schlick's approximation */
                
                R0 = (n1 - n2)/(n1 + n2);
                R0 *= R0;
                double  X = 0.0;
//                if(n1 > n2){
//                    X = 1.0 - sqrtf(1.0 - sin2theta_t);
//                }
//                else{ X = 1.0 - costheta_v; }
                X = 1.0 - costheta_v;
                Refl = R0 + (1.0 - R0) *  X * X * X * X * X;
                Trans = 1.0 - Refl;
                
            }
        }
        else {/* Total internal reflection */
            Refl = 1.0f;
        }
        
        /* Calculate reflection due to reflectance */
        if(bounceCount >0){
            N = hInfo.N;
            Point3 V = Point3(iRay.p.x -  P.x, iRay.p.y - P.y, iRay.p.z - P.z);
            //V.Normalize();
            Point3 VR = 2 * V.Dot(N) * N - V;
            //VR.Normalize();
            Ray rRay = Ray(P, VR);
            //rRay.dir.Normalize();
            rRay.p.x = rRay.p.x + bias *rRay.dir.x;
            rRay.p.y = rRay.p.y + bias *rRay.dir.y;
            rRay.p.z = rRay.p.z + bias *rRay.dir.z;
            HitInfo temp1;
            temp1.Init();
            if(rootNode.GetNumChild()>0){
                if(RayTrace_2(rRay, temp1)){
                    bounceCount --;
                    reflShade =   temp1.node->GetMaterial()->Shade(rRay, temp1, lights, bounceCount);
                }
            }
        }
        
//        cout<<"Refl: "<<Refl<<"Trans "<<Trans<<endl;
        tShade = refraction * (Trans * tShade + Refl * reflShade);
        
        
    }
    





    /* calculate reflection*/
    if(reflection.Grey()>0 && bounceCount > 0){

        Point3 N = hInfo.N;
        Point3 V = Point3(iRay.p.x -  P.x, iRay.p.y - P.y, iRay.p.z - P.z);
       // V.Normalize();
        Point3 VR = 2 * V.Dot(N) * N - V;
        Ray rRay = Ray(hInfo.p, VR);
        //rRay.dir.Normalize();
        rRay.p.x = rRay.p.x + bias *rRay.dir.x;
        rRay.p.y = rRay.p.y + bias *rRay.dir.y;
        rRay.p.z = rRay.p.z + bias *rRay.dir.z;
        HitInfo temp;
        temp.Init();
        if(rootNode.GetNumChild()>0){
            if(RayTrace_2(rRay, temp)){
                bounceCount--;
                rShade = reflection * temp.node->GetMaterial()->Shade(rRay, temp, lights, bounceCount);
            }
        }
    }
    
  
    
    /* Add shade with reflected and refracted colors */
    shade += (rShade + tShade);
    return shade;
};


bool RayTrace_2(const Ray &ray, HitInfo &hitInfo)
{
//    cout<<"Ray trace "<<endl;
    return TraceNode(ray, hitInfo, rootNode);
  
}

bool TraceNode(const Ray &r, HitInfo &hInfo,const Node &node)
{

    Ray ray = r;
    ray = node.ToNodeCoords(r);

    const Object* obj = node.GetObject();
    bool objHitTest = false;
    bool childHit = false;
    
    if(obj)    {
        if(obj->IntersectRay(ray, hInfo)){
            objHitTest=true;
            hInfo.node = &node;
        }
    }
    if (node.GetNumChild() > 0)    {
        for (int i = 0; i < node.GetNumChild(); i++)        {
            Ray r = ray;
            const Node &childNode = *node.GetChild(i);
            if(TraceNode(r, hInfo, childNode)){
                    childHit = true;
            }
        }
        if(childHit)        {
            objHitTest = true;
        }
    }
    
    if(objHitTest){
       node.FromNodeCoords(hInfo);
    }
    return objHitTest;
}

void StopRender()
{
  std::cout<<"Stop Render ..."<<std::endl;
}

void MtlBlinn::SetViewportMaterial() const{
    ColorA c;
    c = diffuse;
    glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r );
    c = specular;
    glMaterialfv( GL_FRONT, GL_SPECULAR, &c.r );
    glMaterialf( GL_FRONT, GL_SHININESS, glossiness );
}

bool RayTrace(HitInfo &hitInfo, Node* curnode, Ray ray, int PixIndex)
{
    Node* node = curnode;
    bool hitTest = false;
    
    const Object *obj = node->GetObject();
    ray = curnode->ToNodeCoords(ray);
    if(obj){
        //       cout<<"Transforming to..."<<endl;
        
        HitInfo tempHitInfo;
        tempHitInfo.Init();
        tempHitInfo.node = node;
        tempHitInfo.z = hitInfo.z;
        hitTest = obj->IntersectRay(ray, tempHitInfo);
        node->FromNodeCoords(tempHitInfo);
        if(hitTest && tempHitInfo.z < hitInfo.z){
            hitInfo = tempHitInfo;
            //cout<<hitInfo.z<<endl;
        }
        //else hitTest=false;
    }
    if(node->GetNumChild()>0)
    {
        //cout<<"Children "<<node->GetNumChild()<<endl;
        for(int i=0;i<curnode->GetNumChild();++i)
        {
            //            cout<<"Child "<<i<<endl;
            node = curnode->GetChild(i);
            HitInfo temp;
            temp.Init();
            temp = hitInfo;
            if(RayTrace(hitInfo, node, ray, PixIndex)){
                curnode->FromNodeCoords(hitInfo);
                
                //                cout<<"Transforming from "<<curnode->GetNumChild()<<endl;
                if(temp.z > hitInfo.z) hitTest = true;
                else{
                    // hitInfo = temp;
                    hitTest = false;
                    continue;
                }
            }
            
        }
    }
    
    
    if(hitTest) return true;
    else return false;
    
}

int main(int argc, char* argv[])
{
  //  int xdim = 400, ydim = 296, zdim = 352; //walnut
  //  int xdim = 512, ydim = 512, zdim = 134; char* datafile="data/Pig.raw";// Pig.raw
  //  int xdim = 512,ydim=512,zdim=63; char* datafile = "data/Teddy_512x512x63.raw";//teddy
  //  int xdim = 256,ydim=256,zdim=161; char* datafile = "data/Tooth_256x256x161.raw";//tooth
  //  int xdim = 256,ydim=256,zdim=512; char* datafile = "data/Carp_256x256x512.raw";//Carp
  //  int xdim = 512, ydim=512,zdim=106; char* datafile = "data/Cadhead_512x512x106.raw";//Cad head
  //  int xdim = 256,ydim=256,zdim=109; char* datafile = "data/MRIwoman_256x256x109.raw";//MRI woman
  //int xdim = 256,ydim=256,zdim=128; char* datafile = "data/engine_256x256x128.raw"; //Engine 8-bit
  //  int xdim = 256,ydim=256,zdim=256; char* datafile = "data/bonsai_256x256x256.raw";
  //  int xdim = 256,ydim=256,zdim=113; char* datafile = "data/CThead_256x256x113.raw";
  //  int xdim = 256,ydim=256,zdim=178; char* datafile = "data/Teapot_256x256x178.raw";
    int xdim = 256,ydim=256,zdim=256; char* datafile = "data/foot_8bit_256x256x256.raw";
  //       int xdim = 128,ydim=256,zdim=256; char* datafile = "data/VisMale_128x256x256.raw";
  //         int xdim = 256,ydim=256,zdim=256; char* datafile = "data/Engine_256x256x256.raw";

  //const char* tf_filename = "data/engine_tf1.1dt";
    const char* tf_filename = "data/foot_2.1dt";
  //const char* tf_filename = "data/VisMale_128x256x256.1dt";
    //  const char* tf_filename = "data/Engine_256x256x256.1dt";
  const char* filename = "scene.xml";
  
  uchar* volumeData = NULL;
  Color* color_tf = NULL;
  float* alpha_tf = NULL;
  VolumeData DataLoader;
  int tf_size = 0;
  uchar minData, maxData;
  LoadScene(filename);
  if(DataLoader.Load(datafile, xdim, ydim, zdim, &volumeData) && DataLoader.LoadTF(tf_filename) ){
    std::cout<<"Loaded data file: "<<datafile<<std::endl;
    theBoxObject.SetDimensions(xdim, ydim, zdim);
    theBoxObject.SetData(volumeData, lights);
    theBoxObject.CalculateGradients();
    DataLoader.GetTransferFunction( &color_tf, &alpha_tf, tf_size, minData, maxData );
    theBoxObject.SetTransferFunction(color_tf, alpha_tf, tf_size, minData, maxData);

    volumeData = NULL;
    color_tf = NULL;
    alpha_tf = NULL;
  }
  else{
    std::cout<<"Failure to load data file"<<std::endl;
    return -1;;
  }

  //glutInit(&argc,argv);
  //ShowViewport();

  return 0;
}

