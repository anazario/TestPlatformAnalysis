#include "PulseShape.h"

PulseShape::PulseShape(vector<vector<float>> sampleCollection){
  samples_ = sampleCollection;
  sampleSize_ = samples_.size();
  stepSize_ = 0.01;
  GetMaxT0Amp0Vec();
}

PulseShape::~PulseShape(){}

void PulseShape::GetMeanErrors(vector<float>& mean, vector<float>& low_err, vector<float>& high_err, int CI){

  int vecSize = histVec_.size();

  float meanVal = -999.;
  float lowErr  = -999.;
  float highErr = -999.;
  
  for(int i = 0; i < vecSize; i++){
    PulseTools::CalcInterval(histVec_[i],CI,meanVal,lowErr,highErr);
    mean.push_back(meanVal);
    low_err.push_back(lowErr);
    high_err.push_back(highErr);
  }
}

void PulseShape::GetErrors(vector<float>& low_err, vector<float>& high_err, float CI){

  int vecSize = histVec_.size();
  float mean;
  float lowErr  = -999.;
  float	highErr	= -999.;

  for(int i = 0; i < vecSize; i++){
    PulseTools::CalcInterval(histVec_[i],CI,mean,lowErr,highErr);
    low_err.push_back(lowErr);
    high_err.push_back(highErr);
  }
}

void PulseShape::GetMaxT0andAmp0(float* sample, float& t0, float& amp0){

  int nSteps = 5000;

  float approx_t0 = float(PulseTools::FindMinAbsolute(sample,SAMPLE_SIZE))/10.;
  float xmin = approx_t0 - 1;
  float xmax = approx_t0 + 1;
  float step_size = (xmax-xmin)/float(nSteps);

  float timetemp = xmin;
  float new_amp = 0;

  t0 = 0;
  amp0 = -999.;

  for(int i = 0; i<nSteps; i++){
    new_amp = -PulseTools::InterpolateFunc(sample,SAMPLE_SIZE,timetemp);
    if(new_amp > amp0){
      amp0 = new_amp;
      t0 = timetemp;
    }
    timetemp += step_size;
  }
}

void PulseShape::GetMaxT0Amp0Vec(){

  float t0, amp0;

  for(int i = 0; i < sampleSize_; i++){
    GetMaxT0andAmp0(&samples_[i][0],t0,amp0);
    tPeak_.push_back(t0);
    maxAmp_.push_back(amp0);
  }
}

float* PulseShape::MakeInterpolation(float* sample, float t0, float amp_peak, float step_size){
    
  int xmin = -1.;
  int xmax = 1.;

  float t0_amp = -PulseTools::InterpolateFunc(sample,SAMPLE_SIZE,t0);
  float LE_ratio = round((t0_amp/amp_peak)*10.)/10.;
  if(LE_ratio != 1.0)
    cout << LE_ratio << endl;
  LE_ratio = 1.;

  //cout << t0 << " " << t0_amp << " " << amp_peak << " " << (t0_amp/amp_peak) << endl;

  //if(LE_ratio != 0.2)
  //cout << "ratio: " << LE_ratio << endl;

  int nLeft, nRight;
  float size = float(xmax-xmin)/step_size;

  if (LE_ratio == 1){
    nLeft = int(size/2.);
    nRight = nLeft;
  }

  else{
    nLeft = int(size*0.20);
    nRight = int(size) - nLeft;
  }

  float arr_size = float(xmax-xmin)/step_size + 1;
  float* pulse = new float[int(arr_size)];

  for(int i = -nLeft; i <= nRight; i++)
    pulse[i+nLeft] = (-PulseTools::InterpolateFunc(sample, SAMPLE_SIZE, t0 + float(i)*step_size)/amp_peak);

  return pulse;
}

float* PulseShape::GetMeanPulse(vector<float> t0, vector<float> max_amp, float step_size, std::function<float*(float*,float,float,float)> PulseType){

  if(histVec_.size()>0)
    histVec_.clear();

  float* processed_sample;

  int vec_size = int(2./step_size)+1;
  histVec_ = PulseTools::MakeHistArr(vec_size,-0.2,1.2,1400);

  for(int i = 0; i < samples_.size(); i++){
    processed_sample = PulseType(&samples_[i][0],t0[i],max_amp[i],step_size);
      
    for(int j = 0; j < vec_size; j++){
      histVec_[j].Fill(processed_sample[j]);
      //if(processed_sample[j] > 1.)
	//cout << processed_sample[j] << endl;
    }
  }
  delete [] processed_sample;
  return PulseTools::GetMeanArr(histVec_);
}

float* PulseShape::GetPulseCDF(float* sample, float t0, float max_amp, float step_size){

  int arr_size = int(2./step_size)+1;
  float* interpolated_sample = MakeInterpolation(sample,t0,max_amp,step_size);
  float* cdf = new float[arr_size];

  float sum = 0;
  float integral = PulseTools::Integral(interpolated_sample,step_size);

  for(int i = 0; i < arr_size; i++){
    if(interpolated_sample[i]>0){
      sum += interpolated_sample[i]*step_size/integral;
    }
    cdf[i] = sum;
  }

  delete [] interpolated_sample;
  return cdf;
}

void PulseShape::PlotHists(){
  if(histVec_.size() > 0){
    TCanvas *cv = new TCanvas("cv","cv",600,800);
    histVec_[0].Draw();
    cv->SaveAs("new_hists/pulse0.pdf");
    for(int i = 1; i < histVec_.size(); i++){
      histVec_[i].Draw();
      cv->SaveAs(Form("new_hists/pulse%i.pdf",i));
    }
  }
  else{
    cout << "Histogram vector is empty!" << endl;
    return;
  }
}

void PulseShape::PlotAllPulses(){
  PlotAll(tPeak_,maxAmp_,stepSize_,"AllPulses","Amplitude (1/Peak Amplitude)",MakeInterpolation);
}

void PulseShape::PlotAllCDFs(){
  PlotAll(tPeak_,maxAmp_,stepSize_,"AllCDFs","Cumulative Distribution",GetPulseCDF);
}

void PulseShape::PlotPulseMean(){
  PlotMean(tPeak_,maxAmp_,stepSize_,"PulseMean","Amplitude (1/Peak Amplitude)",MakeInterpolation);
}

void PulseShape::PlotCDFMean(){
  PlotMean(tPeak_,maxAmp_,stepSize_,"CDFMean","Cumulative Distribution",GetPulseCDF);
}

void PulseShape::PlotPulseMeanErr(){
  PlotMeanErr(tPeak_,maxAmp_,stepSize_,"PulseMeanErr","Amplitude (1/Peak Amplitude)","pulse",MakeInterpolation);
}

void PulseShape::PlotCDFMeanErr(){
  PlotMeanErr(tPeak_,maxAmp_,stepSize_,"CDFMeanErr","Cumulative Distribution","CDF",GetPulseCDF);
}
/*
void PulseShape::PlotRangeCDF(float xmin, float xmax, int total_slices){
  PlotRange(xmin,xmax,total_slices,"rangeplots/CDFMeanErr","Cumulative Distribution","CDF",GetPulseCDF);    
}

void PulseShape::PlotRangePulse(float xmin, float xmax, int total_slices){
  PlotRange(xmin,xmax,total_slices,"rangeplots/PulseMeanErr","Amplitude (1/Peak Amplitude)","pulse",MakeInterpolation);
}
*/
void PulseShape::PlotAll(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, 
				std::function<float*(float*,float,float,float)> PulseType){
  
  int total = int(2./step_size)+1;
  float* xAxis = PulseTools::LinSpace(-1.,1.,step_size);
  float* yAxis;
  yAxis = PulseType(&samples_[0][0],t0[0],max_amp[0],step_size);

  TCanvas *cv = new TCanvas("cv","cv",600,800);
  TGraph *tg = new TGraph[samples_.size()];
  tg[0] = TGraph(total,xAxis,yAxis);
  gStyle->SetOptTitle(kFALSE);
  tg[0].GetXaxis()->SetTitle("time (ns)");
  tg[0].GetYaxis()->SetTitle(y_label);

  tg[0].Draw();

  for(int i = 1; i < sampleSize_; i++){
    yAxis = PulseType(&samples_[i][0],t0[i],max_amp[i],step_size);
    tg[i] = TGraph(total,xAxis,yAxis);
    tg[i].Draw("same");
  }

  CMSmark("");
  gPad->SetGrid(1, 1); gPad->Update();
  cv->SaveAs(name+".pdf");
  cv->Close();

  delete cv;
  delete [] tg;
  delete [] xAxis;
  delete [] yAxis;
}

void PulseShape::PlotMean(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, 
				 std::function<float*(float*,float,float,float)> PulseType){

  int total = int(2./step_size)+1;
  float* xAxis = PulseTools::LinSpace(-1.,1.,step_size);
  float* yAxis = GetMeanPulse(t0,max_amp,step_size,PulseType);

  TCanvas *cv = new TCanvas("cv","cv",600,800);
  TGraph *tg = new TGraph(total,xAxis,yAxis);

  gStyle->SetOptTitle(kFALSE);
  tg->SetLineWidth(2);
  tg->GetXaxis()->SetTitle("time (ns)");
  tg->GetYaxis()->SetTitle(y_label);

  tg->Draw();
  CMSmark("");
  gPad->SetGrid(1, 1); gPad->Update();
  cv->SaveAs(name+".pdf");
  cv->Close();

  delete cv;
  //delete tg;
  delete [] xAxis;
  delete [] yAxis;
}

void PulseShape::PlotMeanErr(vector<float> t0, vector<float> max_amp, float step_size, TString name, TString y_label, TString opt, 
		 std::function<float*(float*,float,float,float)> PulseType){

  int total = int(2./step_size)+1;
  float* xAxis = PulseTools::LinSpace(-1.,1.,step_size);
  float* mean = GetMeanPulse(t0,max_amp,step_size,PulseType);

  vector<float> low_err68;
  vector<float> high_err68;

  vector<float> low_err95;
  vector<float> high_err95;

  vector<float> low_err99;
  vector<float> high_err99;

  GetErrors(low_err68,high_err68,0.68);
  GetErrors(low_err95,high_err95,0.95);
  GetErrors(low_err99,high_err99,0.99);

  TCanvas *cv = new TCanvas("cv","cv",600,800);
  TGraphAsymmErrors* tg = new TGraphAsymmErrors(total,xAxis,mean,0,0,&low_err68[0],&high_err68[0]);
  TGraphAsymmErrors* tg95 = new TGraphAsymmErrors(total,xAxis,mean,0,0,&low_err95[0],&high_err95[0]);
  TGraphAsymmErrors* tg99 = new TGraphAsymmErrors(total,xAxis,mean,0,0,&low_err99[0],&high_err99[0]);
  TGraph* line = new TGraph(total,xAxis,mean);

  //Aesthetics                                                                                                                             
  int style = 3144;
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetPalette(kAvocado);
  gStyle->SetHatchesLineWidth(2);

  line->SetLineWidth(3);
  line->SetLineStyle(1);
  tg->SetFillStyle(style);
  tg95->SetFillStyle(style);
  tg99->SetFillStyle(style);
  tg99->GetXaxis()->SetTitle("time (ns)");
  tg99->GetYaxis()->SetTitle(y_label);
  tg99->GetYaxis()->SetTitleOffset(1.25);
  tg99->SetTitle("99% CL");
  tg95->SetTitle("95% CL");
  tg->SetTitle("68% CL");
  line->SetTitle("Mean Pulse");
  tg99->Draw("a4l PFC");
  tg95->Draw("samel4 PFC");
  tg->Draw("samel4 PFC");
  line->Draw("sameL");

  //Legend
  TLegend* leg;

  if(opt == "pulse")
    leg = new TLegend(0.7,0.75,0.89,0.89);
  else
    leg = new TLegend(0.20,0.7,0.37,0.84);

  leg->AddEntry(line,"Mean Pulse","l");
  leg->AddEntry(tg,"68% CL","f");
  leg->AddEntry(tg95,"95% CL","f");
  leg->AddEntry(tg99,"99% CL","f");
  //leg->SetBorderSize(0);                                                                                                                 
  leg->Draw("same");

  CMSmark("");
  gPad->SetGrid(1, 1); gPad->Update();

  cv->SaveAs(name+".pdf");
  cv->Close();

  delete [] xAxis;
  delete [] mean;
  /*
  delete [] low_err68;
  delete [] high_err68;
  delete [] low_err95;
  delete [] high_err95;
  delete [] low_err99;
  delete [] high_err99;
  */
  delete leg;
  delete tg;
  delete tg95;
  delete tg99;
  delete line;
  delete cv;
}
/*
void PulseShape::PlotRange(float xmin, float xmax, int total_slices, TString name, TString y_label, TString opt, 
				  std::function<float*(float*,float,float,float)> PulseType){

  TString gen_string = "_Amp_%iTo%i";

  vector<float> total_t0;
  vector<float> total_max_amp;
  vector<float> t0;
  vector<float> max_amp;

  //SetCuts(xmin,xmax);
  //SetGoodPulses();

  GetMaxT0Amp0Vec(total_t0,total_max_amp);
  std::sort(total_max_amp.begin(),total_max_amp.end());
  int npulses = total_max_amp.size();

  float amp_percent = (1/float(total_slices));

  float temp_xmin = total_max_amp[0];
  float temp_xmax = total_max_amp[min(npulses-1,int(npulses/total_slices))];

  for(int i = 1; i < total_slices+1; i++){

    SetCuts(temp_xmin,temp_xmax);
    SetGoodPulses();

    cout << "Plotting range " << temp_xmin << " to " << temp_xmax << endl;

    GetMaxT0Amp0Vec(t0,max_amp);

    cout << amp_percent*(i+1) << endl;
    cout << round(npulses*amp_percent)*(i+1) << endl;
 
    //PlotAll(t0,max_amp,step_size,name+Form(gen_string,temp_xmin,temp_xmax),y_label,PulseType);
    //PlotMean(t0,max_amp,step_size,name+Form(gen_string,temp_xmin,temp_xmax),y_label,PulseType);
    PlotMeanErr(t0,max_amp,step_size,name+Form(gen_string,int(round(temp_xmin)),int(round(temp_xmax))),y_label,opt,PulseType);

    temp_xmin = temp_xmax;
    temp_xmax = total_max_amp[min(npulses-1,int(npulses*amp_percent)*(i+1))];

    t0.clear();
    max_amp.clear();
  }
}
*/
