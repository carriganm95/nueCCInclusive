#pragma once

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "sbnana/CAFAna/Core/Var.h"
#include "sbnana/CAFAna/Core/Cut.h"
#include "sbnana/CAFAna/Core/Tree.h"
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

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TVector3.h"
#include "TF1.h"

#include <vector>
#include <variant>

#include "TStyle.h" ///// GM

Color_t color_nue   = kBlue-7;
Color_t color_numu  = kOrange+2;
Color_t color_nc    = kGreen+2;
Color_t color_cos   = kGray+2;
Style_t line_nue    = kSolid;
Style_t line_numu   = kSolid;
Style_t line_nc     = kDashed;
Style_t line_cos    = kSolid;

const double e_mass = 0.00051099895000; // GeV/c^2 //https://pdg.lbl.gov/2025/tables/contents_tables.html

namespace ana {

  extern const Var kTrkLen;
  extern const MultiVar kTrkLenAll;
  extern const Var kLeadingProtonId;
  extern const Var kLeadingMuonId;
  extern const Var kLargestShowerId;
  extern const Var kLongestTrackId;
  extern const Var kLeadingProtonMomentum;
  extern const MultiVar kProtonMomentum;
  extern const Var kLeadingMuonMomentum;
  extern const MultiVar kMuonMomentum;
  extern const Var kTruth_NeutrinoE;
  extern const MultiVar kShowerE;
  extern const Var kLeadingShowerE;

  extern const Cut kShwContainedFD;
  extern const Cut kNotClearCosmic;
  extern const Cut kIsClearCosmic;
  extern const Cut kRFiducialCut;
  extern const Cut kNucrlongtrkdiryCut;
  extern const Cut kLeadingProtonCut;
  extern const Cut kLeadingMuonCut;
  extern const Cut kLargestShowerCut;
  extern const Cut kLongestTrackCut;
  extern const Cut kGoodMuon;
  extern const Cut kGoodProton;
  extern const Cut kGoodEle;

  extern const TruthCut kTFiducialCut;

  float trkDistance(caf::Proxy<caf::SRPFP>&, caf::SRSliceProxy*);
  bool trackFit(caf::Proxy<caf::SRPFP>&);
  bool containedTrack(caf::Proxy<caf::SRPFP>&, float);
  bool containedVertex(const caf::Proxy<caf::SRVector3D>&, float, float, float, float, float, float);
  bool containedShower(caf::Proxy<caf::SRPFP>&);
  bool checkMuon(caf::Proxy<caf::SRPFP>&, const caf::SRSliceProxy*);
  bool checkProton(caf::Proxy<caf::SRPFP>&, const caf::SRSliceProxy*);
  bool checkEle(caf::Proxy<caf::SRPFP>&, const caf::SRSliceProxy*, int);
  bool checkShower(caf::Proxy<caf::SRPFP>&, const caf::SRSliceProxy*);

  static const caf::SRSpillProxy* g_CurrentSpill = nullptr;
  static int g_MaxShowerSliceIdx = -1;

}