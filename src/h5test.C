#include "DataReaderH5.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <fstream>
#include <map>

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
  DataReaderH5 d(("data/Scan15/hdf5/x"+to_string(x)+"_y"+to_string(y)+".hdf5").c_str(),1);
  d.PrintInfo();
  float** samples = d.GetWaveForms();
  float** time = d.GetTimeArr();

  int s = 16;
  
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
  /*
  map<string, int> m;
  
  int value;
  vector<char> keyVec;
  fstream file;
  string  filename, key;
  char c;
  filename = "data/Scan15/Scan15_info.json";
  
  // opening file
  file.open(filename.c_str());
  
  // extracting words from the file
  if(file >> c && c == '{'){
    while (file >> c)
      {
	if (c == '\"'){
	  while (file >> c && c != '\"')
	    key += c;
	}
	if(c == ':' && file >> value){
	  m.insert(pair<string,int>(key,value));
	  key = "";
	}
      }
  }
  cout << m["nTrig"] << endl;
  cout << m["xMaxIdx"] << endl;
  cout << m["yMaxIdx"] << endl;
  */
  /*
  int cnt = 0;
  int idx = 0;
  string str1 = "data/Scan15/hdf5/x0_y0.hdf5";
  string str2 = "data/Scan15";  

  string path;
  string dir;
  
  int l1 = str1.length();
  int l2 = str2.length();

  idx = str1.find_last_of('/');
  if(idx+1 == l1){
    str1.replace(idx,1,"");
    idx = str1.find_last_of('/');
    l1 = str1.length();
  }

  for (int i = 0; i < idx+1; i++)
    path += str1[i];
  for (int i = idx+1; i < l1; i++)
    dir += str1[i];

  cout << path << endl;
  cout << dir << endl;
  
    
  cout << idx << endl;
  idx = str2.find_last_of('/');
  cout << idx << endl;
    */

  return 0;
}

