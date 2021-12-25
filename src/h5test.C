#include "DataReaderH5.h"

int main(void){

  DataReaderH5 d;
  d.GetChData("data/Scan15/hdf5/x18_y16.hdf5",1);
  d.PrintInfo();

  return 0;
}
