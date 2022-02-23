#include <stdio.h>  
#include <unistd.h>

#include "Scan.h"
#include "DataReaderH5.h"
#include "PulseList.h"

#include "TH1.h"
#include "TH2.h"
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
  TH1D *ampHist = new TH1D(("amp_"+name).c_str(),("amp_"+name).c_str(),52,0.,780.);
  TH1D *noiseHist = new TH1D(("noise_"+name).c_str(),("noise_"+name).c_str(),52,0.,780.);
  TH1D *timeHist = new TH1D(("time_"+name).c_str(),("time_"+name).c_str(),44,-120.,320.);
  TH2D *ampVsTime = new TH2D(("ampVsTime_"+name).c_str(),("ampVsTime_"+name).c_str(),44,-120.,320.,52,0.,780.);
  
for(auto file : scan.GetFileList()){

  count++;
  fprintf(stdout, "\r  Processed files: %5d of %5d ", count, totalFiles);
  fflush(stdout);
      
  DataReaderH5 *data = new DataReaderH5((dataLoc+file).c_str());
  double rate = data->GetRate(channel);
  
  if(rate > 0.){
    PulseList pulseList = data->GetPulseList(channel);
    for(int s = 0; s < pulseList.size(); s++){
      noiseHist->Fill(pulseList[s]->GetNoiseRMS()*1e3);
      if(data->GetYpos() > 7. && data->GetYpos() < 20){
	
	double maxAmp = pulseList[s]->GetMaxAmp()*1e3;
	double maxTime = pulseList[s]->GetMaxTime();
	
	ampHist->Fill(maxAmp);
	timeHist->Fill(maxTime);
	ampVsTime->Fill(maxTime,maxAmp);
	
      }
      delete pulseList[s];
    }
  }
  delete data;
 }

 cout << endl;

  Plot1D(ampHist, "plots/amplitude_"+name, "Amplitude [mV]", "Events", true);
  Plot1D(noiseHist, "plots/noise_"+name, "Amplitude [mV]", "Events", true);
  Plot1D(timeHist, "plots/time_"+name, "Time [ns]", "Events", true);
  Plot2D(ampVsTime, "plots/ampVsTime_"+name, "Time [ns]", "Amplitude [mV]",true);

  delete ampHist;
  delete noiseHist;
  delete timeHist;
  delete ampVsTime;
  
  return 0;
}
