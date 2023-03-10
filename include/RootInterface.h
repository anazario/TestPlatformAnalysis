#ifndef ROOTINTERFACE_H
#define ROOTINTERFACE_H

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class RootInterface {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  Double_t        horiz_offset;
  Double_t        horiz_scale;
  Long64_t        samples[1002];
  Double_t        trig_offset;
  Double_t        trig_time;
  Double_t        vert_offset;
  Double_t        vert_scale;

  RootInterface(TTree *tree=0, const int channel=1);
  virtual ~RootInterface();

  int GetSampleSize();
private:

  int sample_size_;

  // List of branches
  TBranch        *b_horiz_offset;
  TBranch        *b_horiz_scale;
  TBranch        *b_samples;
  TBranch        *b_trig_offset;
  TBranch        *b_trig_time;
  TBranch        *b_vert_offset;
  TBranch        *b_vert_scale;
  
  virtual void Init(TTree *tree, const int channel);
  
};

#endif
