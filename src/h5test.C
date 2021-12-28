#include "DataReaderH5.h"
#include "PulseShape.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <fstream>
#include <map>

int main(void){
  /*
  vector<int> v = {1,2,3};
  int *a = new int[3];

  a = &v[0];

  for (int i = 0; i < v.size(); i++)
    cout << a[i] << endl;
  */
  /*
  for(int ch = 1; ch < 5; ch++){
    cout << "Channel: " << ch << endl;
    for(int x = 0; x < 21; x++){
      for(int y = 0; y < 21; y++){
	DataReaderH5 d(("data/Scan15/hdf5/x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),ch);
	if(d.GetRate() != 0 && d.GetRate() > 80.){
	  //cout <<	"xidx: " << x << ", yidx: " << y << endl;
	  d.PrintInfo();
	}
      }
    }
  }
  */
  
  int x = 8;
  int y = 5;
  DataReaderH5 d(("data/Scan15/hdf5/x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),1);
  d.PrintInfo();

  vector<vector<float>> samples = d.GetWaveForms();
  vector<vector<float>> time = d.GetTimeArr();
  
  PulseShape ps(samples);
  ps.PlotAllPulses();
  ps.PlotPulseMeanErr();
  ps.PlotHists();

  /*
  cout << time.size() << endl;
  for(int i = 0; i < 2002; i++)
    cout << time[0][i] << endl;
  cout << samples.size() << endl;
  cout << samples[0].size() << endl;
  */
  /*
  int s = 12;

  vector<float> xVec = time[s];
  vector<float> yVec = samples[s];

  cout << xVec.size() << endl;
  cout << time[s].size() << endl;
  
  TCanvas *cv = new TCanvas("cv","cv",600,800);
  TGraph *g = new TGraph(SAMPLE_SIZE,&xVec[0],&yVec[0]);
  g->SetMarkerStyle(8);
  g->SetMarkerSize(0.5);
  g->Draw("AP");
  cv->SaveAs("graph.pdf");
  */

  return 0;
}

