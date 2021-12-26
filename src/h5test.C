#include "DataReaderH5.h"
#include "TGraph.h"
#include "TCanvas.h"

int main(void){
  /*
  for(int ch = 1; ch < 5; ch++){
    cout << "Channel: " << ch << endl;
    for(int x = 0; x < 21; x++){
      for(int y = 0; y < 21; y++){
	DataReaderH5 d(("data/Scan15/hdf5/x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),ch);
	if(d.GetRate() != 0 && d.GetRate() > 50.){
	  cout <<	"xidx: " << x << ", yidx: " << y << endl;
	  d.PrintInfo();
	}
      }
    }
  }
  */
  
  int x = 10;
  int y = 5;
  DataReaderH5 d(("data/Scan15/hdf5/x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),4);
  d.PrintInfo();
  float** samples = d.GetWaveForms();
  float** time = d.GetTimeArr();//PulseTools::LinSpace(0.,1.0,1./SAMPLE_SIZE);

  //cout << time[500] << endl;

  int s = 11;
  
  TCanvas *cv = new TCanvas("cv","cv",600,800);
  TGraph *g = new TGraph(SAMPLE_SIZE,time[s],samples[s]);
  g->SetMarkerStyle(8);
  g->SetMarkerSize(0.5);
  g->Draw("AP");
  cv->SaveAs("graph.pdf");
  /*
  for(int i = 0; i < N_TRIG; i++){
    for(int j = 0; j < SAMPLE_SIZE; j++){
      cout << samples[i][j] << endl;
    }
  }
  */
  return 0;
}
