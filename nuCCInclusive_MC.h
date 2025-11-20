#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Cuts/TruthCuts.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/SRVector3D.h"

#include "sbnana/SBNAna/Vars/Vars.h"
#include "sbnana/SBNAna/Vars/Binnings.h"
#include "sbnana/SBNAna/Vars/NumuVars.h"
//#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Vars/TruthVars.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "sbnana/SBNAna/Cuts/TruthCuts.h"
//#include "sbnana/SBNAna/Vars/NumuVarsIcarus202106.h"

#include "sbnana/SBNAna/Vars/NueVars.h"
#include "sbnana/SBNAna/Cuts/NueCuts.h"
#include "nuCCInclusiveCuts.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEfficiency.h"

#include <vector>
#include <variant>

#include "TStyle.h" ///// GM

using namespace ana;

using Var      = ana::_Var<caf::Proxy<caf::SRSlice>>;
using MultiVar = ana::_MultiVar<caf::Proxy<caf::SRSlice>>;

// Note: monostate is first, and IS default-constructible
using VarOrMulti = std::variant<std::monostate, Var, MultiVar>;

struct HistAxisVM {
  // X axis
  std::string labelX;
  ana::Binning binsX;
  VarOrMulti   varX;

  // Y axis (optional; if is2D == false, these are dummy)
  std::string labelY;
  ana::Binning binsY;
  VarOrMulti   varY;

  bool is2D;  // flag to tell your code whether to treat it as 1D or 2D

  // 1D constructor
  HistAxisVM(const std::string &lab,
             const ana::Binning &b,
             const VarOrMulti &v)
    : labelX(lab),
      binsX(b),
      varX(v),
      labelY(""),
      // use a *valid*, but trivial, binning for Y (e.g. 1 bin) – you must pass it in
      binsY(ana::Binning::Simple(1, 0., 1.)), // or another simple public ctor
      varY(std::monostate{}),
      is2D(false)
  {}

  // 2D constructor
  HistAxisVM(const std::string &lab1,
             const ana::Binning &b1,
             const VarOrMulti &v1,
             const std::string &lab2,
             const ana::Binning &b2,
             const VarOrMulti &v2)
    : labelX(lab1),
      binsX(b1),
      varX(v1),
      labelY(lab2),
      binsY(b2),
      varY(v2),
      is2D(true)
  {}
};

template <class C>
class cutFlow {
  public:
      // Data members
      std::string name_;
      std::vector<C> cuts_; //!                   
      std::vector<std::string> cutNames_;
      
      std::vector<C> nm1Cuts_; //!                   
      std::vector<TH1*> nm1H_;                     
      std::vector<std::unique_ptr<Spectrum>> nm1S_;
      std::vector<std::string> nm1Names_;
      std::vector<std::string> nm1PltVars_;

      std::vector<C> seqCuts_; //!                   
      std::vector<TH1*> seqH_;                     
      std::vector<std::unique_ptr<Spectrum>> seqS_;
      std::vector<std::string> seqNames_;
      std::vector<std::string> seqPltVars_;

      std::vector<TH1*> hFinal_;

      double POT = 3.0E20;
      
      std::vector<Color_t> colors = {kBlack, kRed, kGreen, kMagenta, kCyan, kGray, kPink+10, kSpring+8, kOrange+7, kGreen+4, kYellow, kBlue, kViolet+7, kOrange+3};

      // Default constructor
      cutFlow() = default;

      // Constructor with only name
      explicit cutFlow(std::string name)
          : name_(std::move(name)) {}

      // Full constructor
      cutFlow(std::string name, std::vector<C> cuts, std::vector<std::string> cutNames)
          : name_(std::move(name)),
            cuts_(std::move(cuts)), 
            cutNames_(std::move(cutNames)) {}


      // Begin member functions

      // Create NM1 cuts from cuts_
      void createNM1Cuts(C initialCut=kNoCut){
        for (size_t t = 0; t < cuts_.size(); ++t) {
          C combined = initialCut;
  
          for (size_t s = 0; s < cuts_.size(); ++s) {
            if (s == t) continue; // skip the t-th cut
            combined = combined && cuts_[s];
          }
          nm1Cuts_.push_back(combined);
        }
      }

      // Create sequential cuts from cuts_
      void createSeqCuts(C initialCut=kNoCut){
        for (size_t t = 0; t < cuts_.size(); ++t) {
          C combined = initialCut;
  
          for (size_t s = 0; s < t+1; ++s) {
            //cout << "t/s: " << t << ", " << s << std::endl;
            combined = combined && cuts_[s];
          }
          seqCuts_.push_back(combined);
        }
      }

      // Create NM1 Spectra objects
      void createNM1Spectra(SpectrumLoader &loader, HistAxisVM ax, std::string name="", SpillCut spillCut=kNoSpillCut){
        nm1PltVars_ = {name};
        for (size_t c = 0; c < nm1Cuts_.size(); ++c){
          std::visit(
            [&](auto const &v) {
              using T = std::decay_t<decltype(v)>;
              if constexpr (std::is_same_v<T, std::monostate>) {
                return;  // or assert, or throw, depending on your policy
              } 
              else {
                nm1S_.emplace_back(
                  std::make_unique<ana::Spectrum>(
                    name, ax.binsX, loader, v, spillCut, nm1Cuts_[c]
                  )
                );
              }
            },
            ax.varX
          );
          TString n = name != "" ? cutNames_[c]+"_"+name : cutNames_[c];
          nm1Names_.emplace_back(n);
        }
      }

      // Create NM1 Spectra objects for multiple variables
      void createNM1Spectra(SpectrumLoader &loader, std::vector<HistAxisVM> ax, std::vector<std::string> name, SpillCut spillCut=kNoSpillCut){
        if (name.size() != ax.size()) cout << "Error: number of variables and names does not match" << endl;
        nm1PltVars_ = name;
        for (size_t i = 0; i < ax.size(); ++i){
          for (size_t c = 0; c < nm1Cuts_.size(); ++c){
            if (ax[i].is2D){
                // 2D: visit both X and Y variants simultaneously
                std::visit(
                [&](auto const &vx, auto const &vy) {
                  using TX = std::decay_t<decltype(vx)>;
                  using TY = std::decay_t<decltype(vy)>;

                  // Skip if either axis is monostate
                  if constexpr (std::is_same_v<TX, std::monostate> || std::is_same_v<TY, std::monostate>) {
                    return;
                  }
                  // Require both axes to be the same concrete variant type (e.g. Var vs Var, MultiVar vs MultiVar)
                  else if constexpr (!std::is_same_v<TX, TY>) {
                    return;
                  }
                  else if constexpr (std::is_same_v<TX, MultiVar> || std::is_same_v<TY, MultiVar>) {
                    return; // skip 2D spectrum creation if either axis is MultiVar
                  }
                  else {
                    // Assume Spectrum has a 2D constructor: (name, binsX, binsY, loader, varX, varY, spillCut, cut)
                    nm1S_.emplace_back(
                      std::make_unique<ana::Spectrum>(
                      name[i], loader, ax[i].binsX, vx, ax[i].binsY, vy, spillCut, nm1Cuts_[c]
                      )
                    );
                  }
                },
                ax[i].varX, ax[i].varY
                );
            }
            else{
              std::visit(
                [&](auto const &v) {
                  using T = std::decay_t<decltype(v)>;
                  if constexpr (std::is_same_v<T, std::monostate>) {
                    // This is the "empty" / default state that ROOT may construct.
                    // You must *not* call Spectrum with it.
                    return;  // or assert, or throw, depending on your policy
                  } else {
                    nm1S_.emplace_back(
                      std::make_unique<ana::Spectrum>(
                        name[i], ax[i].binsX, loader, v, spillCut, nm1Cuts_[c]
                      )
                    );
                }
                },
                ax[i].varX
              );
            } 
            nm1Names_.emplace_back(TString(cutNames_[c]+"_"+name[i]));
          }
        }
      }

      // Create NM1 Spectra objects for multiple truth variables
      void createNM1Spectra(SpectrumLoader &loader, std::vector<std::string> labels, std::vector<Binning> bins, std::vector<TruthVar> vars, std::vector<std::string> name, TruthCut truthCut, SpillCut spillCut=kNoSpillCut){
        if (name.size() != labels.size()) cout << "Error: number of variables and names does not match" << endl;
        nm1PltVars_ = name;
        for (size_t i = 0; i < labels.size(); ++i){
          for (size_t c = 0; c < nm1Cuts_.size(); ++c){
            nm1S_.emplace_back(std::make_unique<Spectrum>(labels[i], bins[i], loader, vars[i], truthCut, spillCut, nm1Cuts_[c]));
            nm1Names_.emplace_back(TString(cutNames_[c]+"_"+name[i]));
          }
        }
      }

      // Create sequential Spectra objects
      void createSeqSpectra(SpectrumLoader &loader, HistAxisVM ax, std::string name="", SpillCut spillCut=kNoSpillCut){
        seqPltVars_ = {name};
        for (size_t c = 0; c < seqCuts_.size(); ++c){
          std::visit(
            [&](auto const &v) {
              using T = std::decay_t<decltype(v)>;
              if constexpr (std::is_same_v<T, std::monostate>) {
                // This is the "empty" / default state that ROOT may construct.
                // You must *not* call Spectrum with it.
                return;  // or assert, or throw, depending on your policy
              } 
              else {
                seqS_.emplace_back(
                  std::make_unique<ana::Spectrum>(
                    name, ax.binsX, loader, v, spillCut, seqCuts_[c]
                  )
                );
              }
            },
            ax.varX
          );
          //seqS_.emplace_back(std::make_unique<Spectrum>(name, ax.bins[0], loader, ax.vm[0], spillCut, seqCuts_[c]));
          //seqS_.emplace_back(std::make_unique<Spectrum>(loader, ax, kNoSpillCut, seqCuts_[c]));
          TString n = name != "" ? cutNames_[c]+"_"+name : cutNames_[c];
          seqNames_.emplace_back(n);
        }
      }

      // Create sequential Spectra objects for multiple variables
      void createSeqSpectra(SpectrumLoader &loader, std::vector<HistAxisVM> ax, std::vector<std::string> name, SpillCut spillCut=kNoSpillCut){
        if (name.size() != ax.size()) cout << "Error: number of variables and names does not match" << endl;
        seqPltVars_ = name;
        for (size_t i = 0; i < ax.size(); ++i){

          for (size_t c = 0; c < seqCuts_.size(); ++c){
            if (ax[i].is2D){
                // 2D: visit both X and Y variants simultaneously
                std::visit(
                [&](auto const &vx, auto const &vy) {
                  using TX = std::decay_t<decltype(vx)>;
                  using TY = std::decay_t<decltype(vy)>;

                  // Skip if either axis is monostate
                  if constexpr (std::is_same_v<TX, std::monostate> || std::is_same_v<TY, std::monostate>) {
                    return;
                  }
                  // Require both axes to be the same concrete variant type (e.g. Var vs Var, MultiVar vs MultiVar)
                  else if constexpr (!std::is_same_v<TX, TY>) {
                    return;
                  }
                  else if constexpr (std::is_same_v<TX, MultiVar> || std::is_same_v<TY, MultiVar>) {
                    return; // skip 2D spectrum creation if either axis is MultiVar
                  }
                  else {
                    // Assume Spectrum has a 2D constructor: (name, binsX, binsY, loader, varX, varY, spillCut, cut)
                    seqS_.emplace_back(
                      std::make_unique<ana::Spectrum>(
                      name[i], loader, ax[i].binsX, vx, ax[i].binsY, vy, spillCut, seqCuts_[c]
                      )
                    );
                  }
                },
                ax[i].varX, ax[i].varY
                );
            }
            else{
              std::visit(
                [&](auto const &v) {
                  using T = std::decay_t<decltype(v)>;
                  if constexpr (std::is_same_v<T, std::monostate>) {
                    // This is the "empty" / default state that ROOT may construct.
                    // You must *not* call Spectrum with it.
                    return;  // or assert, or throw, depending on your policy
                  } else {
                    seqS_.emplace_back(
                      std::make_unique<ana::Spectrum>(
                        name[i], ax[i].binsX, loader, v, spillCut, seqCuts_[c]
                      )
                    );
                }
                },
                ax[i].varX
              );
            }
            seqNames_.emplace_back(TString(cutNames_[c]+"_"+name[i]));
          }
        }
      }

      // Create NM1 Spectra objects for multiple truth variables
      void createSeqSpectra(SpectrumLoader &loader, std::vector<std::string> labels, std::vector<Binning> bins, std::vector<TruthVar> vars, std::vector<std::string> name, TruthCut truthCut, SpillCut spillCut=kNoSpillCut){
        if (name.size() != labels.size()) cout << "Error: number of variables and names does not match" << endl;
        seqPltVars_ = name;
        for (size_t i = 0; i < labels.size(); ++i){
          for (size_t c = 0; c < seqCuts_.size(); ++c){
            seqS_.emplace_back(std::make_unique<Spectrum>(labels[i], bins[i], loader, vars[i], truthCut, spillCut, seqCuts_[c]));
            seqNames_.emplace_back(TString(cutNames_[c]+"_"+name[i]));
          }
        }
      }

      // Create TH1s for all nm1 spectra
      void createNM1Hist(){
        for (auto& s : nm1S_){
          nm1H_.emplace_back(s->ToTHX(POT));
        }
      }

      // Create TH1s for all sequential spectra
      void createSeqHist(){
        for (auto& s : seqS_){
          seqH_.emplace_back(s->ToTHX(POT));
        }
        //hFinal_ = (TH1*)seqH_.back()->Clone(TString("hFinal_"+name_));
      }

      // Save NM1 Plots
      void saveNM1Plots(std::string plotDir="nm1Plots"){
        gDirectory->mkdir(plotDir.c_str());
        gDirectory->cd(plotDir.c_str());
        for (size_t h=0; h < nm1H_.size(); h++) {
          nm1H_[h]->Write(TString("nm1"+nm1Names_[h]));
        }
        gDirectory->cd("../");
      }

      // Save NM1 Plots
      void saveSeqPlots(std::string plotDir="seqPlots"){
        gDirectory->mkdir(plotDir.c_str());
        gDirectory->cd(plotDir.c_str());
        for (size_t h=0; h < seqH_.size(); h++) {
          seqH_[h]->Write(TString("seq"+seqNames_[h]));
        }
        gDirectory->cd("../");
      }

      void createCombinedNM1(){
        TCanvas* c = new TCanvas("c", "c");
        float maxY = 0;
        TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9);
        for (size_t i=0; i < nm1Names_.size(); ++i){
          TString drawS = "hist";
          if (i > 0) drawS += " same";

          //get the max/min y of the hists
          if(i % nm1Cuts_.size() == 0){
            for (size_t j=i; j < i+nm1Cuts_.size(); ++j){
              if (nm1H_[j]->GetMaximum() > maxY) maxY = nm1H_[j]->GetMaximum();
            }
            maxY = maxY*1.1;
            nm1H_[i]->GetYaxis()->SetRangeUser(0, maxY);
          }

          // change hist colors and draw
          nm1H_[i]->SetLineColor(colors[i%nm1Cuts_.size()]);
          nm1H_[i]->Draw(drawS);

          // add cut to legend
          if (i/nm1Cuts_.size() == 0) l->AddEntry(nm1H_[i], TString(cutNames_[i]), "l");

          // write out the canvas after cuts are all added for a given var
          if (i > 0 && ((i+1) % nm1Cuts_.size()) == 0){
            l->Draw();
            c->Write(TString("nm1"+nm1PltVars_[i/nm1Cuts_.size()]));
            c->Clear();
            maxY = 0;
          }
        }
      }

      void createCombinedSeq(){
        TCanvas* c = new TCanvas("c", "c");
        float maxY = 0;
        TLegend *l = new TLegend(0.7, 0.6, 0.9, 0.9);
        for (size_t i=0; i < seqNames_.size(); ++i){
          TString drawS = "hist";
          if (i > 0) drawS += " same";

          //get the max/min y of the hists
          if(i % seqCuts_.size() == 0){
            for (size_t j=i; j < i+seqCuts_.size(); ++j){
              if (seqH_[j]->GetMaximum() > maxY) maxY = seqH_[j]->GetMaximum();
            }
            maxY = maxY*1.1;
            seqH_[i]->GetYaxis()->SetRangeUser(0, maxY);
          }

          // change hist colors and draw
          seqH_[i]->SetLineColor(colors[i%seqCuts_.size()]);
          seqH_[i]->Draw(drawS);

          // add cut to legend
          if (i/seqCuts_.size() == 0) l->AddEntry(seqH_[i], TString(cutNames_[i]), "l");

          // write out the canvas after cuts are all added for a given var
          if (i > 0 && ((i+1) % seqCuts_.size()) == 0){
            l->Draw();
            c->Write(TString("seq"+seqPltVars_[i/seqCuts_.size()]));
            c->Clear();
            maxY = 0;
          }
        }
      }

      std::vector<std::reference_wrapper<TH1*>> getFinalHists() {
          
          std::vector<std::reference_wrapper<TH1*>> result;

          int numHists = seqH_.size() / seqCuts_.size();
          result.reserve(numHists);

          std::vector<int> indices = {};
          for (size_t i=0; i < seqH_.size(); ++i){
            if (i % seqCuts_.size() == 0) indices.push_back(i);
          }

          for (int idx : indices) {
              if (idx < 0 || idx >= static_cast<int>(seqH_.size())) {
                  throw std::out_of_range("Index out of range");
              }
              result.push_back(seqH_[idx]); // stores reference
          }

          std::cout << "size of final hists " << result.size() << endl; 

          return result;
      }

};

template <class C>
void compareCutflows(TFile& f, std::vector<cutFlow<C>*>& cutFlows, std::vector<std::string>& names, std::string dir=""){

  //check to make sure all cutflows have the same number of variables plotted
  int vars = 0;
  int cuts = 0;
  for(size_t i=0; i < cutFlows.size(); ++i){
    if (i == 0){
      vars = cutFlows[i]->seqPltVars_.size();
      cuts = cutFlows[i]->seqCuts_.size();
    }
    if ((int)cutFlows[i]->seqPltVars_.size() != vars || (int)cutFlows[i]->seqCuts_.size() != cuts){
      cout << "The number of variables plotted and cuts applied does not match between the cutflows" << endl;
      return;
    }
  }

  //save plot of each var to the same canvas
  TCanvas* c = new TCanvas("c", "c");
  TLegend* l = new TLegend(0.7, 0.7, 0.9, 0.9);
  int iC = -1;
  for(int v=0; v < vars; ++v){
    c->Clear();
    iC = (v*cuts)+cuts-1;

    float maxY = 0;
    for(size_t i=0; i < cutFlows.size(); ++i){
      if (cutFlows[i]->seqH_[iC]->GetMaximum() > maxY) maxY = cutFlows[i]->seqH_[iC]->GetMaximum();
    }
    maxY *= 1.1;

    for(size_t i=0; i < cutFlows.size(); ++i){
      TString drawS = "hist";
      if (i == 0) cutFlows[i]->seqH_[iC]->GetYaxis()->SetRangeUser(0, maxY);
      if (i > 0) drawS += " same";
      cutFlows[i]->seqH_[iC]->SetLineColor(cutFlows[i]->colors[i]);
      cutFlows[i]->seqH_[iC]->Draw(drawS);
      if (v == 0) l->AddEntry(cutFlows[i]->seqH_[iC], TString(names[i]), "l");
    }
    l->Draw();
    f.cd();
    if (dir != "") gDirectory->cd(dir.c_str());
    c->Write(TString(cutFlows[0]->seqPltVars_[v]));
  }

  return;
}

template <class C>
void compareCutflowsNM1(TFile& f, std::vector<cutFlow<C>*>& cutFlows, std::vector<std::string>& names, std::string cutName, std::string cutVar, std::string dir=""){

  //check to make sure all cutflows have the same number of variables plotted
  int vars = 0;
  int cuts = 0;
  for(size_t i=0; i < cutFlows.size(); ++i){
    if (i == 0){
      vars = cutFlows[i]->nm1PltVars_.size();
      cuts = cutFlows[i]->nm1Cuts_.size();
    }
    if ((int)cutFlows[i]->nm1PltVars_.size() != vars || (int)cutFlows[i]->nm1Cuts_.size() != cuts){
      cout << "The number of variables plotted and cuts applied does not match between the cutflows" << endl;
      return;
    }
  }

  auto cutFind = std::find(cutFlows[0]->cutNames_.begin(), cutFlows[0]->cutNames_.end(), cutName);
  auto varFind = std::find(cutFlows[0]->nm1PltVars_.begin(), cutFlows[0]->nm1PltVars_.end(), cutVar);

  int cutId = cutFind != cutFlows[0]->cutNames_.end() ? std::distance(cutFlows[0]->cutNames_.begin(), cutFind) : -1;
  int varId = varFind != cutFlows[0]->nm1PltVars_.end() ? std::distance(cutFlows[0]->nm1PltVars_.begin(), varFind) : -1;
  if (cutId < 0 || varId < 0) {
    cout << "could not find the cut and/or the variable in compareCutflowNM1" << endl;
    cout << "cutId " << cutId << ", varId " << varId << endl;
    return;
  }

  //save plot of each var to the same canvas
  TCanvas* c = new TCanvas("c", "c");
  TLegend* l = new TLegend(0.7, 0.7, 0.9, 0.9);
  int iC = (varId*cuts)+cutId;

  float maxY = 0;
  for(size_t i=0; i < cutFlows.size(); ++i){
    if (cutFlows[i]->nm1H_[iC]->GetMaximum() > maxY) maxY = cutFlows[i]->nm1H_[iC]->GetMaximum();
  }
  maxY = maxY * 1.1;

  for(size_t i=0; i < cutFlows.size(); ++i){
    TString drawS = "hist";
    if (i > 0) drawS += " same";
    else cutFlows[i]->nm1H_[iC]->GetYaxis()->SetRangeUser(0, maxY);
    cutFlows[i]->nm1H_[iC]->SetLineColor(cutFlows[i]->colors[i]);
    cutFlows[i]->nm1H_[iC]->Draw(drawS);
    l->AddEntry(cutFlows[i]->nm1H_[iC], TString(names[i]), "l");
  }
  l->Draw();
  f.cd();
  if (dir != "") gDirectory->cd(dir.c_str());
  c->Write(TString("nm1" + cutFlows[0]->nm1PltVars_[varId]));

  return;
}



