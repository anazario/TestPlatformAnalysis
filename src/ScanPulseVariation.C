#include <stdio.h>  
#include <unistd.h>

#include "Scan.h"
#include "DataReaderH5.h"
#include "PulseList.h"
#include "PulseVariation.h"

#include "TH1D.h"
#include "TCanvas.h"
#include "Plot.h"

using namespace std;

int main(int argc, char *argv[]){

  int opt;
  int channel = -1;

  char* pathname;
  
  while((opt = getopt(argc, argv, "p:c:")) != -1){
    switch(opt)  
      {
      case 'p':
	pathname = optarg;
	break;
      case 'c':
	channel = atoi(optarg);
	break;	
      }
  }

  string path = string(pathname);

  if(path.find_last_of('/') == path.length()-1)
    path = (path+"hdf5/").c_str();
  else
    path = (path+"/hdf5/").c_str();

  if(!std::__fs::filesystem::is_directory(path)){
    cout << "'" << path << "' is not a scan results directory!" << endl;
    return 0;
  }
  
  if(channel < -1 || channel > 4){
    cout << "Channel must be between an integer between 1 and 4!" << endl;
    return 0;
  }

  /*
  DataReaderH5 dataReader((path+"x0_y0.hdf5").c_str(),1);
  
  int xMaxIdx = dataReader.GetXmaxIdx() + 1;
  int yMaxIdx = dataReader.GetYmaxIdx() + 1;
  int totalFiles = xMaxIdx*yMaxIdx;
  int count = 1;
  */

  Scan scan(path);

  int totalFiles = scan.GetTotalFiles();
  int count = 0;

  string dataLoc = scan.GetPathToData();
  
  cout << "Processing channel: " << channel << endl;
    
  string name = ("ch"+to_string(channel)+"_hist").c_str();

  int interpolationSize = 200;
  double start = -10.;
  double end = 30.;
  vector<TH1D> pulseVariationHists;

  vector<double> interpolationTime = LinSpaceVec(start,end,interpolationSize);
  
  for(int h = 0; h < interpolationSize; h++){
    TH1D hist(TString((name+"_"+to_string(h)).c_str()),TString((name+"_"+to_string(h)).c_str()),14000,-0.2,1.2);
    pulseVariationHists.push_back(hist);
  }
  
  for(auto file : scan.GetFileList()){
  
    fprintf(stdout, "\r  Processed files: %5d of %5d ", count, totalFiles);
    fflush(stdout);
    
    //DataReaderH5 *data = new DataReaderH5((path+"x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),channel);
    DataReaderH5 *data = new DataReaderH5((dataLoc+file).c_str());
    double rate = data->GetRate(channel);
    
    if(rate > 0.){
      PulseList pulseList = data->GetPulseList(channel);
      
      if(data->GetYpos() > 7. && data->GetYpos() < 20){
	for(int s = 0; s < pulseList.size(); s++){
	  Pulse *pulse = pulseList[s];
	  double maxAmp = pulse->GetMaxAmp()*1e3;
	  double maxTime = pulse->GetMaxTime();
	  
	  if(maxAmp > 100. && maxTime < 120.){
	    //vector<double> interpolatedPulse = pulse->GetInterpolation(interpolationSize,start,end,maxTime);
	    vector<double> interpolatedPulse = pulse->GetCDF(interpolationSize,start,end,maxTime);
	    
	    for(int h = 0; h < interpolatedPulse.size(); h++){
	      pulseVariationHists[h].Fill(interpolatedPulse[h]);
	    }
	  }
	  delete pulse;
	}
      }
    }
    delete data;
    count++;

  }
  cout << endl;

  PulseVariation pulseVar(pulseVariationHists,interpolationTime);
  //pulseVar.PlotMeanErrors("PulseMeanErr_Scan15_Ch2","Amplitude/(Max Amplitude)");
  pulseVar.PlotMeanErrors("CDFMeanErr_ScanCh2_1","Pulse Cumulative Distribution");
  pulseVar.PlotHistograms();

  return 0;
}
