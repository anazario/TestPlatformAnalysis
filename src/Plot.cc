#include "Plot.h"

using namespace std;

void CMSmark(TString plotTitle){
  TLatex l;
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextSize(0.035);
  l.SetTextFont(42);
  l.DrawLatex(0.51,0.91,plotTitle);
  l.SetTextSize(0.04);
  l.SetTextFont(42);
  l.DrawLatex(0.11,0.91,"#bf{CMS} #it{Preliminary}");
}

void Plot1D(TH1D* hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog){

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(11111111);
  TCanvas can("can_"+name,"can_"+name,600,800);
  if(isLog)
    can.SetLogy();

  can.SetLeftMargin(0.12);
  can.SetRightMargin(0.3);
  can.SetBottomMargin(0.1);
  can.Update();
  
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(2);
  hist->Draw();

  can.SaveAs(name+".pdf");
  can.Close();
}

void Plot2D(TH2D *hist, const TString name, const TString xlabel, const TString ylabel, const bool isLog){

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11111111);
  
  TCanvas can("can_"+name,"can_"+name,600,800);
  if(isLog)
    can.SetLogz();

  can.SetLeftMargin(0.12);
  can.SetRightMargin(0.2);
  can.SetBottomMargin(0.1);
  can.Update();
  hist->GetXaxis()->SetTitle(xlabel);
  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetZaxis()->SetTitle("Events");
  hist->GetZaxis()->SetTitleOffset(1.5);

  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->CenterTitle();
  hist->GetZaxis()->CenterTitle();
  
  hist->Draw("colz");
  can.SaveAs(name+".pdf");
  can.Close();
}

void PlotGraph(TString name, TString xlabel, TString ylabel, int size, std::vector<double> xVec, std::vector<double> yVec){

  TGraph tg(size,&xVec[0],&yVec[0]);

  TCanvas can("can_"+name,"can_"+name,600,800);
  tg.SetLineWidth(2);
  tg.GetXaxis()->SetTitle(xlabel);
  tg.GetYaxis()->SetTitle(ylabel);
  tg.Draw("AL");
  tg.SetTitle("");
  can.SetLeftMargin(0.15);
  can.SaveAs(name+".pdf");

  can.Close();
}

void PlotWaveform(const std::vector<double>& time, const std::vector<double>& waveform, const int waveform_id) {
    // Create a TGraph object from the input data
    TGraph* graph = new TGraph(time.size(), &time[0], &waveform[0]);

    // Set the axis labels and titles
    graph->GetXaxis()->SetTitle("time [ns]");
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->SetTitle("Amplitude [V]");
    graph->GetYaxis()->CenterTitle();

    // Create a canvas to draw the plot
    TCanvas* canvas = new TCanvas("canvas", "canvas", 800, 600);
    canvas->SetGrid();

    // Draw the TGraph object on the canvas as a line plot
    graph->Draw("AL");

    // Save the plot to a file
    canvas->SaveAs(("waveform"+to_string(waveform_id)+".pdf").c_str());

    // Clean up memory
    delete graph;
    delete canvas;
}


TGraph* PlotScatter(TString name, TString xlabel, TString ylabel, std::vector<double> xVec, std::vector<double> yVec){

  int size = xVec.size();
  
  TGraph* tg = new TGraph(size,&xVec[0],&yVec[0]);

  //tg->SetMarkerStyle(8);
  tg->GetXaxis()->SetTitle(xlabel);
  tg->GetYaxis()->SetTitle(ylabel);
  tg->SetMarkerSize(0.8);
  //tg->Draw("AP");

  return tg;
}

void makeScatterPlot(const std::vector<double> &x_signal, const std::vector<double> &x_background,
                     const std::vector<double> &y_signal, const std::vector<double> &y_background,
                     std::string name, std::string xlabel, std::string ylabel) {
    // Create a new TCanvas and set its title
    TCanvas* canvas = new TCanvas(name.c_str(), name.c_str(), 800, 600);
    canvas->SetTitle(name.c_str());

    // Create a TGraph object for the signal data and set its color to green
    TGraph* signal = new TGraph(x_signal.size(), &x_signal[0], &y_signal[0]);
    signal->SetMarkerStyle(20);
    signal->SetMarkerSize(0.8);
    signal->SetMarkerColor(kGreen);

    // Create a TGraph object for the background data and set its color to blue
    TGraph* background = new TGraph(x_background.size(), &x_background[0], &y_background[0]);
    background->SetMarkerStyle(20);
    background->SetMarkerSize(0.8);
    background->SetMarkerColor(kBlue);

    // Create a legend and add entries for the signal and background TGraphs
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetFillColor(kWhite);
    legend->SetBorderSize(0);
    legend->AddEntry(signal, "Signal", "p");
    legend->AddEntry(background, "Background", "p");

    // Draw the signal and background TGraphs on the canvas
    signal->Draw("AP");
    background->Draw("P same");

    // Set the labels for the x and y axes
    signal->GetXaxis()->SetTitle(xlabel.c_str());
    signal->GetYaxis()->SetTitle(ylabel.c_str());

    // Add the legend to the canvas
    legend->Draw();

    // Update the canvas and save it as a PDF file
    canvas->Update();
    canvas->SaveAs((name + ".pdf").c_str());

    // Free memory allocated for TGraphs, TCanvas, and TLegend objects
    delete signal;
    delete background;
    delete legend;
    delete canvas;
}
