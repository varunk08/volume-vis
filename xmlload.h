#ifndef XMLLOAD_H
#define XMLLOAD_H

#include "scene.h"
#include "objects.h"
#include "sphere.h"
#include "TriObj.h"
#include "Plane.h"
#include "materials.h"
#include "lights.h"
#include "tinyxml/tinyxml.h"
#include "boxobject.h"

int LoadScene(const char *filename);
void LoadScene(TiXmlElement *element);
void LoadNode(Node *node, TiXmlElement *element, int level=0);
void LoadTransform( Transformation *trans, TiXmlElement *element, int level );
void LoadMaterial(TiXmlElement *element);
void LoadLight(TiXmlElement *element);
void ReadVector(TiXmlElement *element, Point3 &v);
void ReadColor (TiXmlElement *element, Color  &c);
void ReadFloat (TiXmlElement *element, float  &f, const char *name="value");


#endif
