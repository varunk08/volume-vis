#ifndef VOLUMEDATA_H
#define VOLUMEDATA_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

bool LoadVolumeData(const char* filename, int xdim, int ydim, int zdim, unsigned short** data);

#endif
