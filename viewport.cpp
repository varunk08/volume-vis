
//-------------------------------------------------------------------------------
///
/// \file       viewport.cpp
/// \author     Cem Yuksel (www.cemyuksel.com)
/// \version    5.0
/// \date       September 23, 2013
///
/// \brief Example source for CS 6620 - University of Utah.
///
//-------------------------------------------------------------------------------
#include "viewport.h"
/*
#include <GLUT/glut.h>
#include <time.h>
#include "scene.h"
#include "objects.h"
*/
//-------------------------------------------------------------------------------
//void Sphere::ViewportDisplay() const
//{
//  static GLUquadric *q = NULL;
//  if ( q == NULL ) q = gluNewQuadric();
//  gluSphere(q,1,50,50);
//}
//void TriObj::ViewportDisplay() const
//{
//  glBegin(GL_TRIANGLES);
//  for ( unsigned int i=0; i<NF(); i++ ) {
//      for ( int j=0; j<3; j++ ) {
//          if ( HasTextureVertices() ) glTexCoord3fv( &VT( FT(i).v[j] ).x );
//          if ( HasNormals() ) glNormal3fv( &VN( FN(i).v[j] ).x );
//          glVertex3fv( &V( F(i).v[j] ).x );
//      }
//  }
//  glEnd();
//}
//void MtlBlinn::SetViewportMaterial() const
//{
//  ColorA c;
//  c = diffuse;
//  glMaterialfv( GL_FRONT, GL_AMBIENT_AND_DIFFUSE, &c.r );
//  c = specular;
//  glMaterialfv( GL_FRONT, GL_SPECULAR, &c.r );
//  glMaterialf( GL_FRONT, GL_SHININESS, shininess );
//}
//void GenLight::SetViewportParam(int lightID, ColorA ambient, ColorA intensity, Point4 pos ) const
//{
//  glEnable ( GL_LIGHT0 + lightID );
//  glLightfv( GL_LIGHT0 + lightID, GL_AMBIENT,  &ambient.r );
//  glLightfv( GL_LIGHT0 + lightID, GL_DIFFUSE,  &intensity.r );
//  glLightfv( GL_LIGHT0 + lightID, GL_SPECULAR, &intensity.r );
//  glLightfv( GL_LIGHT0 + lightID, GL_POSITION, &pos.x );
//}
//-------------------------------------------------------------------------------

void BeginRender(); // Called to start rendering (renderer must run in a separate thread)
void StopRender();  // Called to end rendering (if it is not already finished)

extern Node rootNode;
extern Camera camera;
extern RenderImage renderImage;
extern LightList lights;

//-------------------------------------------------------------------------------

enum Mode {
    MODE_READY,         // Ready to render
    MODE_RENDERING,     // Rendering the image
    MODE_RENDER_DONE    // Rendering is finished
};

enum ViewMode
{
    VIEWMODE_OPENGL,
    VIEWMODE_IMAGE,
    VIEWMODE_Z
};

static Mode     mode        = MODE_READY;       // Rendering mode
static ViewMode viewMode    = VIEWMODE_OPENGL;  // Display mode
static int      startTime;                      // Start time of rendering

//-------------------------------------------------------------------------------

void GlutDisplay();
void GlutReshape(int w, int h);
void GlutIdle();

void GlutMouse(int button, int state, int x, int y);

//-------------------------------------------------------------------------------

void ShowViewport()
{
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );
    if (glutGet(GLUT_SCREEN_WIDTH) > 0 && glutGet(GLUT_SCREEN_HEIGHT) > 0){
        glutInitWindowPosition( (glutGet(GLUT_SCREEN_WIDTH) - camera.imgWidth)/2, (glutGet(GLUT_SCREEN_HEIGHT) - camera.imgHeight)/2 );
    }
    else glutInitWindowPosition( 50, 50 );
    glutInitWindowSize(camera.imgWidth, camera.imgHeight);
    
    glutCreateWindow("Ray Tracer - Volume Visualization");
    glutDisplayFunc(GlutDisplay);
    glutReshapeFunc(GlutReshape);
    glutIdleFunc(GlutIdle);
    glutKeyboardFunc(GlutKeyboard);
    
    glClearColor(0,0,0,0);
    
    glPointSize(3.0);
    glEnable( GL_CULL_FACE );
    
    float zero[] = {0,0,0,0};
    glLightModelfv( GL_LIGHT_MODEL_AMBIENT, zero );
    
    glEnable(GL_NORMALIZE);
    
    glLineWidth(2);
    //hack to start the application
    GlutKeyboard(' ', 0, 0);
    glutMainLoop();
}

//-------------------------------------------------------------------------------

void GlutReshape(int w, int h)
{
    if( w != camera.imgWidth || h != camera.imgHeight ) {
        glutReshapeWindow( camera.imgWidth, camera.imgHeight);
    } else {
        glViewport( 0, 0, w, h );
        
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        float r = (float) w / float (h);
        gluPerspective( camera.fov, r, 0.02, 1000.0);
        
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
    }
}

//-------------------------------------------------------------------------------

void DrawNode( Node *node )
{
    glPushMatrix();
    
    const Material *mtl = node->GetMaterial();
    if ( mtl ) mtl->SetViewportMaterial();
    
    Matrix3 tm = node->GetTransform();
    Point3 p = node->GetPosition();
    float m[16] = { tm[0],tm[1],tm[2],0, tm[3],tm[4],tm[5],0, tm[6],tm[7],tm[8],0, p.x,p.y,p.z,1 };
    glMultMatrixf( m );
    
    Object *obj = node->GetObject();
    if ( obj ) obj->ViewportDisplay();
    
    for ( int i=0; i<node->GetNumChild(); i++ ) {
        DrawNode( node->GetChild(i) );
    }
    
    glPopMatrix();
}

//-------------------------------------------------------------------------------

void DrawScene(bool capture=false)
{

    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
    
    glEnable( GL_LIGHTING );
    glEnable( GL_DEPTH_TEST );
    
    glPushMatrix();
    Point3 p = camera.pos;
    Point3 t = camera.pos + camera.dir;
    Point3 u = camera.up;
    gluLookAt( p.x, p.y, p.z,  t.x, t.y, t.z,  u.x, u.y, u.z );
    
    for ( unsigned int i=0; i<lights.size(); i++ ) {
        lights[i]->SetViewportLight(i);
    }
    
    DrawNode(&rootNode);
    
    glPopMatrix();
    
    glDisable( GL_DEPTH_TEST );
    glDisable( GL_LIGHTING );
    
    if ( capture ) {
       glReadPixels( 0, 0, camera.imgWidth, camera.imgHeight, GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels() );
    }
}

//-------------------------------------------------------------------------------

void DrawProgressBar()
{
    int rp = renderImage.GetNumRenderedPixels();
    int np = renderImage.GetWidth() * renderImage.GetHeight();
    if ( rp >= np ) return;
    
    float done = (float) rp / (float) np;
    
    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();
    
    glBegin(GL_LINES);
    glColor3f(1,1,1);
    glVertex2f(-1,-1);
    glVertex2f(done*2-1,-1);
    glColor3f(0,0,0);
    glVertex2f(done*2-1,-1);
    glVertex2f(1,-1);
    glEnd();
    
    glPopMatrix();
    glMatrixMode( GL_MODELVIEW );
}

//-------------------------------------------------------------------------------

void GlutDisplay()
{
    switch ( viewMode ) {
        case VIEWMODE_OPENGL:
            DrawScene();
            break;
        case VIEWMODE_IMAGE:
            glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
	    glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_RGB, GL_UNSIGNED_BYTE, renderImage.GetPixels() );
            DrawProgressBar();
            break;
        case VIEWMODE_Z:
            glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT );
            if ( ! renderImage.GetZBufferImage() ) renderImage.ComputeZBufferImage();
            glDrawPixels( renderImage.GetWidth(), renderImage.GetHeight(), GL_LUMINANCE, GL_UNSIGNED_BYTE, renderImage.GetZBufferImage() );
            break;
    }
    
    glutSwapBuffers();
}

//-------------------------------------------------------------------------------

void GlutIdle()
{
    static int lastRenderedPixels = 0;
    if ( mode == MODE_RENDERING ) {
        int nrp = renderImage.GetNumRenderedPixels();
	printf("Nrp: %d\n",nrp);
        if ( lastRenderedPixels != nrp ) {
            lastRenderedPixels = nrp;
            if ( renderImage.IsRenderDone() ) {
	      printf("Render done");
                mode = MODE_RENDER_DONE;
                int endTime = (int) time(NULL);
                int t = endTime - startTime;
                int h = t / 3600;
                int m = (t % 3600) / 60;
                int s = t % 60;
                printf("\nRender time is %d:%02d:%02d.\n",h,m,s);
                
                printf("Render Complete");
	    }
            glutPostRedisplay();
        }
    }
}

//-------------------------------------------------------------------------------

void GlutKeyboard(unsigned char key, int x, int y)
{
    switch ( key ) {
        case 27:    // ESC
            exit(0);
            break;
        case ' ':
            switch ( mode ) {
                case MODE_READY:
                    mode = MODE_RENDERING;
                    viewMode = VIEWMODE_IMAGE;
                    DrawScene(true);
                    startTime = time(NULL);
                    BeginRender();
                    break;
                case MODE_RENDERING:
                    mode = MODE_READY;
                    StopRender();
                    glutPostRedisplay();
                    break;
                case MODE_RENDER_DONE:
                    mode = MODE_READY;
                    viewMode = VIEWMODE_OPENGL;
                    glutPostRedisplay();
                    break;
            }
            break;
        case '1':
            viewMode = VIEWMODE_OPENGL;
            glutPostRedisplay();
            break;
        case '2':
            viewMode = VIEWMODE_IMAGE;
            glutPostRedisplay();
            break;
        case '3':
            viewMode = VIEWMODE_Z;
            glutPostRedisplay();
        break;  }
}

//-------------------------------------------------------------------------------
