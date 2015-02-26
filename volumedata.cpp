#include "volumedata.h"

bool LoadVolumeData(const char* filename, int xdim, int ydim, int zdim, unsigned short** data)
{

  FILE *pFile;
  const int size = xdim*ydim*zdim;
  
  pFile = fopen(filename, "rb");
  if(NULL == pFile)
    {
      return false;
    }
  unsigned short *pVolume = new unsigned short[size];
  fread(pVolume,sizeof(unsigned short),size,pFile);
  fclose(pFile);
  /*int index = 0;
    for(int z=0; z < zdim; z++)
      {
	for(int y = 0; y < ydim; y++)
	  {
	    for(int x = 0; x < xdim; x++)
	      {
		//		DataPoints[x][y][z] =(int) theData[index];
		std::cout<< pVolume[index]<<std::endl;
		index++;
	      }
	  }
	
      }
  */

 *data = pVolume;
  return true;

}
