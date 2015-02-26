#include "scene.h"
#include <pthread.h>
pthread_mutex_t mutex1=PTHREAD_MUTEX_INITIALIZER;

void RenderImage::PutPixel(int index, Color24 color,float zbuf)
    {
        img[index]=color;
        zbuffer[index] = zbuf;
        //std::cout<<"zbuf: "<<index<<std::endl;
        pthread_mutex_lock(&mutex1);
        IncrementNumRenderPixel(1);
        pthread_mutex_unlock(&mutex1);
    }
