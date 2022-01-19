#include <stdio.h>  
#include <unistd.h>

#include "DataReaderH5.h"
#include "Pulse.h"
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

  DataReaderH5 dataReader((path+"x0_y0.hdf5").c_str(),1);
  
  int xMaxIdx = dataReader.GetXmaxIdx() + 1;
  int yMaxIdx = dataReader.GetYmaxIdx() + 1;
  int totalFiles = xMaxIdx*yMaxIdx;
  int count = 1;
  
  cout << "Processing channel: " << channel << endl;
    
  string name = ("ch"+to_string(channel)+"_hist").c_str();

  int interpolationSize = 500;
  vector<TH1F> pulseVariationHists;

  vector<float> interpolationTime;
  
  for(int h = 0; h < interpolationSize; h++){
    TH1F hist(TString((name+"_"+to_string(h)).c_str()),TString((name+"_"+to_string(h)).c_str()),14000,-0.2,1.2);
    pulseVariationHists.push_back(hist);
  }
  
  for(int x = 0; x < xMaxIdx; x++){
    for(int y = 0; y < yMaxIdx; y++){

      fprintf(stdout, "\r  Processed files: %5d of %5d ", count, totalFiles);
      fflush(stdout);
      
      DataReaderH5 *dR = new DataReaderH5((path+"x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),channel);
      float rate = dR->GetRate();

      if(rate > 0.){
	vector<vector<float>> samples = dR->GetWaveForms();
	vector<vector<float>> time = dR->GetTimeArr();
	if(dR->GetYpos() > 7. && dR->GetYpos() < 20){
	  for(int s = 0; s < samples.size(); s++){
	    Pulse *pulse = new Pulse(samples[s],time[s]);
	    float maxAmp = pulse->GetMaxAmp()*1e3;
	    float maxTime = pulse->GetMaxTime();
	    
	    if(maxAmp > 100. && maxTime < 120.){
	      vector<float> interpolatedPulse;
	      
	      //pulse->GetInterpolatedPulse(interpolationSize,interpolatedPulse,interpolationTime);
	      pulse->GetCDFpulse(interpolationSize,interpolatedPulse,interpolationTime);
	      
	      for(int h = 0; h < interpolatedPulse.size(); h++){
		pulseVariationHists[h].Fill(interpolatedPulse[h]);
	      }
	    }
	    delete pulse;
	  }
	}
      }
      delete dR;
      count++;
    }
  }
  cout << endl;

  PulseVariation pulseVar(pulseVariationHists,interpolationTime);
  //pulseVar.PlotMeanErrors("CDFMeanErr_Scan15_Ch2","Amplitude/(Max Amplitude)");
  pulseVar.PlotMeanErrors("CDFMeanErr_Scan15_Ch2","Pulse Cumulative Distribution");
  pulseVar.PlotHistograms();

  return 0;
}
