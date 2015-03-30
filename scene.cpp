#include "scene.h"

void RenderImage::PutPixel(int index, Color24 color,float zbuf)
    {
        img[index]=color;
        zbuffer[index] = zbuf;
        IncrementNumRenderPixel(1);

    }
