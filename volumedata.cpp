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
  CreateHistogram(pVolume, size);
 *data = pVolume;
  return true;

}


void CreateHistogram(unsigned short* volume_data, int size)
{
  std::cout<<"Creating histogram..."<<std::endl;
  int min = 999999;
  int max = 0;
  for(int i=0; i<size; i++){
    if((int)volume_data[i] > max) { max=(int)volume_data[i]; }
    if((int)volume_data[i] < min) { min=(int)volume_data[i]; }
  }
  std::cout<<"Min: "<<min<<std::endl;
  std::cout<<"Max: "<<max<<std::endl;
  int binVal = 1000;
  int n = max / binVal;
  int* bins = new int[n+1]();
  
  for( int i=0; i<size; i++ ){
    for( int j=0; j<=n; j++){
      if( volume_data[i] >= j * binVal && volume_data[i] < (j+1) * binVal){
	bins[j]++;
	continue;
      }
      else if( j == n && volume_data[i] > j * binVal){
	bins[j]++;
      }
    }
  }
int sum = 0;
 int max_bin = 0;
 int max_i = -1;
  for(int i=0; i<=n; i++){
    sum += bins[i];
    if( bins[i] > max_bin ){
      max_bin = bins[i];
      max_i = i;
    }
    //std::cout<<"Bin "<<i<<": "<<bins[i]<<std::endl;
  }
  //std::cout<<"SUM: "<<sum<<std::endl;
  //std::cout<<"Max bin index: "<<max_i<<std::endl;
  std::cout<<"Max bin value: "<<max_bin<<std::endl;
  delete[] bins;
  return;
}

