#include "RootInterface.h"

int get_integer_from_string(const std::string& str){
    // Find the first and last position of square brackets
    size_t start_pos = str.find('[');
    size_t end_pos = str.find(']');

    if (start_pos == std::string::npos || end_pos == std::string::npos || start_pos >= end_pos) {
        // No square brackets found or they are in the wrong order
        // Return an error value, like -1
        return -1;
    }

    // Extract the integer value from the string inside the square brackets
    std::stringstream ss(str.substr(start_pos + 1, end_pos - start_pos - 1));
    int integer_value;
    ss >> integer_value;

    if (ss.fail() || !ss.eof()) {
        // Failed to parse an integer from the string inside the square brackets
        // Return an error value, like -1
        return -1;
    }

    // Return the parsed integer value
    return integer_value;
}


std::vector<std::string> get_branch_names(TTree* tree) {
    std::vector<std::string> branch_names;
    TObjArray* branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* branch = static_cast<TBranch*>(branches->At(i));
        branch_names.push_back(std::string(branch->GetName()));
    }
    return branch_names;
}

RootInterface::RootInterface(TTree *tree, const int channel) : fChain(0){
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("signal_test_ch1.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("signal_test_ch1.root");
    }
    f->GetObject("data",tree);    
  }

  for(auto &str : get_branch_names(tree)){
    sample_size_ = get_integer_from_string(str);
    if(sample_size_ > 0)
      break;
  }

  Init(tree, channel);
}

RootInterface::~RootInterface(){
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

void RootInterface::Init(TTree *tree, const int channel){
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress(("ch"+std::to_string(channel)+"_horiz_offset").c_str(), &horiz_offset, &b_horiz_offset);
  fChain->SetBranchAddress(("ch"+std::to_string(channel)+"_horiz_scale").c_str(), &horiz_scale, &b_horiz_scale);
  fChain->SetBranchAddress(("ch"+std::to_string(channel)+"_samples["+std::to_string(sample_size_)+"]").c_str(), samples, &b_samples);
  fChain->SetBranchAddress(("ch"+std::to_string(channel)+"_trig_offset").c_str(), &trig_offset, &b_trig_offset);
  fChain->SetBranchAddress(("ch"+std::to_string(channel)+"_trig_time").c_str(), &trig_time, &b_trig_time);
  fChain->SetBranchAddress(("ch"+std::to_string(channel)+"_vert_offset").c_str(), &vert_offset, &b_vert_offset);
  fChain->SetBranchAddress(("ch"+std::to_string(channel)+"_vert_scale").c_str(), &vert_scale, &b_vert_scale);
}

int RootInterface::GetSampleSize(){
  return sample_size_;
}

