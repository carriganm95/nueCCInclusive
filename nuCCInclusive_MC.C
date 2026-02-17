#pragma once 

#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Core/Spectrum.h"

#include "nuCCInclusiveCuts.h"
#include "nuCCInclusiveCuts.C"
#include "nuCCInclusive_MC.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TEfficiency.h"

using namespace ana;

void nuCCInclusive_MC(bool debug=false, std::string jobId="", std::string outputFile="debug.root", int start=-1, int stop=-1){

  std::cout << "jobId: " << jobId << " , start: " << start << ", stop: " << stop << std::endl;
  TString fphase2 = "/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/00/*/*.flat.caf*.root";
  if ( debug ) fphase2 = "/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/05/04/detsim_2d_icarus_refactored_detsim_stage0_stage1_79315342_0.flat.caf-0985b172-de17-487e-8183-915eed082119.root";
  else if (jobId != "") fphase2 = jobId;

  std::cout << "Using input: " << fphase2 << std::endl;

  // option to use CVN scores
  bool useCVN = false;
  
  // option to use nugraph scores
  bool useNuGraph = false;

  if (useNuGraph){
    SetTrackScoreBranch("ngscore");
    std::cout << "Using track score branch: " << trackScoreBranch << std::endl;
  }

  std::string inputStr = fphase2.Data();

  // Check if input is a file list (.txt or .list)
  bool isFileList = (inputStr.size() > 4 && 
                    (inputStr.substr(inputStr.size()-4) == ".txt" || 
                      inputStr.substr(inputStr.size()-5) == ".list"));

  SpectrumLoader loader = [&]() {
    if (isFileList) {
      // Read file list and build vector of files
      std::vector<std::string> allFiles;
      std::ifstream inFile(inputStr);
      if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file list: " << inputStr << std::endl;
        exit(1);
      }
      std::string line;
      while (std::getline(inFile, line)) {
        // Skip empty lines and comments
        if (!line.empty() && line[0] != '#') {
          allFiles.push_back(line);
        }
      }
      inFile.close();
      
      // Extract subset based on start/stop indices
      std::vector<std::string> fileList;
      int startIdx = (start >= 0) ? start : 0;
      int stopIdx = (stop > 0) ? stop : allFiles.size();
      
      // Validate indices
      if (startIdx >= (int)allFiles.size()) {
        std::cerr << "Error: start index " << startIdx << " exceeds file list size " << allFiles.size() << std::endl;
        exit(1);
      }
      stopIdx = std::min(stopIdx, (int)allFiles.size());
      
      fileList.assign(allFiles.begin() + startIdx, allFiles.begin() + stopIdx);
      std::cout << "Loaded " << fileList.size() << " files (indices " << startIdx 
                << " to " << stopIdx << ") from list: " << inputStr << std::endl;
      return SpectrumLoader(fileList);
    }
    else {
      // Check for wildcards
      bool hasWildcards = (inputStr.find("?") != std::string::npos || 
                          inputStr.find("*") != std::string::npos);
      return hasWildcards 
        ? SpectrumLoader(inputStr)
        : SpectrumLoader(std::vector<std::string>{inputStr});
    }
  }();

  
  double POT = 3.0E20; //BNB 6.6 NuMI 6
  const double nominalIntensity =  5.13E13; // POT per spill, 5.13 one data (run 8515) when we get a stable beam 5E13 //tested with 6E13
  const double nominalLivetime = POT/nominalIntensity;

  TFile f(TString(outputFile), "RECREATE");

  const Cut totalCut = kNoCut;

  std::vector<Cut> finalCuts;
  std::vector<std::string> finalCutNames;

  // Add fiducial and barycenter cuts to all cutflows
  finalCuts.insert(finalCuts.end(), {totalCut, kLargestShowerCut, kNotClearCosmic, kRFiducialCut, kBarycenterCuts}); 
  finalCutNames.insert(finalCutNames.end(), {"NoCut", "Shower", "ClearCosmic", "Fiducial", "Barycenter"});

  // if using CVN cut
  if ( useCVN ) {
    finalCuts.push_back(kCVNNueCut);
    finalCutNames.push_back("CVNNue");
  }

  // track length and shower containment cuts
  finalCuts.insert(finalCuts.end(), {kLongestTrackLCut, kShwContainedCut});
  finalCutNames.insert(finalCutNames.end(), {"LongestTrack", "ShwContained"});

  // if using shower energy selection
  //finalCuts.push_back(kShowerEnergy);
  //finalCutNames.push_back("ShowerEnergy");

  // shower quality cuts
  finalCuts.insert(finalCuts.end(), {kShowerDedx, kShowerConversionGap, kShowerDensity});
  finalCutNames.insert(finalCutNames.end(), {"ShowerDedx", "ShowerConversionGap", "ShowerDensity"});

  std::vector<Cut> eShowerCuts = finalCuts;
  std::vector<std::string> cutNames = finalCutNames;
  
  //std::vector<Cut> eShowerCuts = {totalCut, kRFiducialCut, kBarycenterCuts, kLongestTrackLCut, kShwContainedCut, kShowerEnergy, kShowerDedx, kShowerConversionGap, kShowerDensity};
  //std::vector<std::string> cutNames = {"NoCut", "Fiducial", "Barycenter", "LongestTrack", "ShwContained", "ShowerEnergy", "ShowerDedx", "ShowerConversionGap", "ShowerDensity"};

  // version with cvn cut, shower energy removed
  //std::vector<Cut> eShowerCuts = {totalCut, kRFiducialCut, kBarycenterCuts, kCVNNueCut, kLongestTrackLCut, kShwContainedCut, kShowerDedx, kShowerConversionGap, kShowerDensity};
  //std::vector<std::string> cutNames = {"NoCut", "Fiducial", "Barycenter", "CVNNue", "LongestTrack", "ShwContained", "ShowerDedx", "ShowerConversionGap", "ShowerDensity"};

  // Add the max shower energy slice cut to the beginning of the cuts for efficiency calculations to ensure only one slice per spill is considered
  auto eShowerCutEff = eShowerCuts;
  auto cutNamesEff = cutNames;
  // eShowerCutEff.insert(eShowerCutEff.begin(), kIsMaxShowerEnergySlice());
  // cutNamesEff.insert(cutNamesEff.begin(), "MaxShowerEnergySlice");
  eShowerCutEff.insert(eShowerCutEff.begin(), kSliceNeutrinoMatched);
  cutNamesEff.insert(cutNamesEff.begin(), "NeutrinoMatchedSlice");

  cutFlow<Cut> all("all", eShowerCuts, cutNames); // all events
  cutFlow<Cut> nue("nue", eShowerCuts, cutNames); // nue CC events
  cutFlow<Cut> numu("numu", eShowerCuts, cutNames); // numu CC events
  cutFlow<Cut> nc("nc", eShowerCuts, cutNames); // NC events
  cutFlow<Cut> cosmic("cosmic", eShowerCuts, cutNames); // cosmic events
  cutFlow<Cut> fid("fid", eShowerCuts, cutNames); // fiducial events
  cutFlow<Cut> nueEff("nueEff", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueP("nueP", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (neutrino only)
  cutFlow<Cut> nueA("nueA", eShowerCutEff, cutNamesEff);  // nue CC events selecting only slice with max shower energy per spill (antineutrino only)
  cutFlow<Cut> nueCCQE("nueCCQE", eShowerCutEff, cutNamesEff); // nue CCQE events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueCCRes("nueCCRes", eShowerCutEff, cutNamesEff); // nue CCRes events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueCCDis("nueCCDis", eShowerCutEff, cutNamesEff);  // nue CCDis events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueCCCoh("nueCCCoh", eShowerCutEff, cutNamesEff); // nue CCCoh events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueCCMEC("nueCCMEC", eShowerCutEff, cutNamesEff); // nue CCMEC events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueNC("nueNC", eShowerCutEff, cutNamesEff); // nue NC events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueCCOther("nueCCOther", eShowerCutEff, cutNamesEff); // nue CCOther events selecting only slice with max shower energy per spill
  cutFlow<Cut> nueEle("nueEle", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (true electron only)
  cutFlow<Cut> nueMu("nueMu", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (true muon only)
  cutFlow<Cut> nuePro("nuePro", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (true proton only)
  cutFlow<Cut> nuePi0("nuePi0", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (true neutral pion only)
  cutFlow<Cut> nuePiC("nuePiC", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (true charged pion only)
  cutFlow<Cut> nuePh("nuePh", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (true photon only)
  cutFlow<Cut> nueOther("nueOther", eShowerCutEff, cutNamesEff); // nue CC events selecting only slice with max shower energy per spill (true other particle only)

  //Cut initialCut = kLargestShowerCut && kNotClearCosmic;
  Cut initialCut = kNoCut;

  all.createNM1Cuts(initialCut);
  all.createSeqCuts(initialCut);

  nue.createNM1Cuts(initialCut && kNueCC);
  nue.createSeqCuts(initialCut && kNueCC);

  numu.createNM1Cuts(initialCut && kNumuCC);
  numu.createSeqCuts(initialCut && kNumuCC);

  nc.createNM1Cuts(initialCut && kNC);
  nc.createSeqCuts(initialCut && kNC);

  cosmic.createNM1Cuts(initialCut && kThisCosmic);
  cosmic.createSeqCuts(initialCut && kThisCosmic);

  fid.createNM1Cuts(initialCut);
  fid.createSeqCuts(initialCut);

  nueEff.createNM1Cuts(initialCut && kNueCC);
  nueEff.createSeqCuts(initialCut && kNueCC);

  nueP.createNM1Cuts(initialCut && kTrueNue);
  nueP.createSeqCuts(initialCut && kTrueNue);

  nueA.createNM1Cuts(initialCut && kTrueANue);
  nueA.createSeqCuts(initialCut && kTrueANue);

  nueCCQE.createNM1Cuts(initialCut && kNueCC && kQECut);
  nueCCQE.createSeqCuts(initialCut && kNueCC && kQECut);

  nueCCRes.createNM1Cuts(initialCut && kNueCC && kResCut);
  nueCCRes.createSeqCuts(initialCut && kNueCC && kResCut);

  nueCCDis.createNM1Cuts(initialCut && kNueCC && kDisCut);
  nueCCDis.createSeqCuts(initialCut && kNueCC && kDisCut);

  nueCCCoh.createNM1Cuts(initialCut && kNueCC && kCohCut);
  nueCCCoh.createSeqCuts(initialCut && kNueCC && kCohCut);

  nueCCMEC.createNM1Cuts(initialCut && kNueCC && kMECCut);
  nueCCMEC.createSeqCuts(initialCut && kNueCC && kMECCut);

  nueNC.createNM1Cuts(initialCut && kIsNue && kNC);
  nueNC.createSeqCuts(initialCut && kIsNue && kNC);

  nueCCOther.createNM1Cuts(initialCut && kNueCC && kOtherModeCut);
  nueCCOther.createSeqCuts(initialCut && kNueCC && kOtherModeCut);

  nueMu.createNM1Cuts(initialCut && kNueCC && kIsTrueMuon);
  nueMu.createSeqCuts(initialCut && kNueCC && kIsTrueMuon);

  nuePro.createNM1Cuts(initialCut && kNueCC && kIsTrueProton);
  nuePro.createSeqCuts(initialCut && kNueCC && kIsTrueProton);

  nuePi0.createNM1Cuts(initialCut && kNueCC && kIsTrueNeutralPion);
  nuePi0.createSeqCuts(initialCut && kNueCC && kIsTrueNeutralPion);

  nuePiC.createNM1Cuts(initialCut && kNueCC && kIsTrueChargedPion);
  nuePiC.createSeqCuts(initialCut && kNueCC && kIsTrueChargedPion);

  nuePh.createNM1Cuts(initialCut && kNueCC && kIsTruePhoton);
  nuePh.createSeqCuts(initialCut && kNueCC && kIsTruePhoton);

  nueOther.createNM1Cuts(initialCut && kNueCC && kOtherTruthPDG);
  nueOther.createSeqCuts(initialCut && kNueCC && kOtherTruthPDG);

  // Spectrum
  ////######################################################/////

  std::vector<std::string> varNames = {"tNueEnergy", "rShowerEnergy", "rShowerDeDx", "rShowerDensity", "rShowerConversionGap", "rLongestTrackLength", 
                                      "rBarycenterDeltaZ", "rBarycenterRadius", "rFiducial", "rShowerContained", "rMuonChi2", "tInteractionMode", "tNeutrinoType",
                                      "rSubleadingShowerEnergy", "rShowerEResidual", "rSecondaryShowerEResidual", "rNeutrinoScore", "tShowerLength", "rShowerLenResidual", "rVertexResidual",
                                      "tShowerE", "tSubleadingShowerE", "kLongestTrackLengthResidual", "tTLongestTrackLength", "kShowerCosAngle", "kShowerP",
                                      "tFiducial", "rShowerLength", "rShowerEndX", "rShowerEndY", "rShowerEndZ", "rShowerStartX", "rShowerStartY", "rShowerStartZ",
                                      "rVertexCryo", "rTruthPDGID", "rClearCosmic", "rLargestShowerCut", "tGoodNue"//, "rAllShowerEnergies"
                                      };

  std::vector<HistAxisVM> ax = {HistAxisVM("True E_{#nu} (GeV)", Binning::Simple(60,0.f,6.f),  kTruth_NeutrinoE),
                              HistAxisVM("Reco Shower E (GeV)", Binning::Simple(60,0.f,6.f),  kLeadingShowerE),
                              HistAxisVM("Reco Shower DeDx (MeV/cm)", Binning::Simple(50, 0.f, 10.f), kLeadingShowerDedx),
                              HistAxisVM("Reco Shower Density (MeV/cm)", Binning::Simple(150, 0.f, 30.f), kLeadingShowerDensity),
                              HistAxisVM("Reco Shower Conversion Gap (cm)", Binning::Simple(50, 0.f, 10.f), kLeadingShowerConversionGap),
                              HistAxisVM("Longest Track Length (cm)", Binning::Simple(40, 0, 400), kLongestTrack),
                              HistAxisVM("Barycenter Delta(Z) (cm)", Binning::Simple(48, 0, 120), kBaryDeltaZT),
                              HistAxisVM("Barycenter Radius (cm)", Binning::Simple(48, 0, 120), kBaryRadiusT),
                              HistAxisVM("Vertex in Fiducial Volume", Binning::Simple(2, 0, 2), kRFiducial),
                              HistAxisVM("Shower Contained in Fiducial Volume", Binning::Simple(2, 0, 2), kShwContained),
                              HistAxisVM("Muon Chi2 Score", Binning::Simple(30, 0, 120), kMuonChi2),
                              HistAxisVM("Genie Mode", Binning::Simple(15, -1, 14), kInteractionMode),
                              HistAxisVM("Neutrino Type", Binning::Simple(15, -1, 14), kNeutrinoType),
                              HistAxisVM("Subleading Shower Energy (GeV)", Binning::Simple(60,0.f,6.f),  kSubLeadingShowerE),
                              HistAxisVM("Shower Energy Residual (GeV)", Binning::Simple(40, -2, 2), kShowerEResidual),
                              HistAxisVM("Shower Energy Residual (GeV)", Binning::Simple(40, -2, 2), kSubleadingShowerEResidual),
                              HistAxisVM("Neutrino Score", Binning::Simple(50, 0, 1), kNuScore),
                              HistAxisVM("True Shower Length (cm)", Binning::Simple(40, 0, 400), kTShowerLength),
                              HistAxisVM("Shower Length Residual (Reco-True)/Reco", Binning::Simple(40, -2, 2), kShowerLenResidual),
                              HistAxisVM("Vertex Residual (Reco-True)/Reco", Binning::Simple(40, -2, 2), kVertexResidual),
                              HistAxisVM("True Leading Shower Energy", Binning::Simple(60, 0, 6), kTShowerE),
                              HistAxisVM("True Subleading Shower Energy", Binning::Simple(60, 0, 6), kTSubleadingShowerE),
                              HistAxisVM("Leading Track Length Residual (Reco-True)/Reco", Binning::Simple(40, -2, 2), kTrackLenResidual),
                              HistAxisVM("True Leading Track Length", Binning::Simple(40, 0, 400), kTTrackLen),
                              HistAxisVM("Leading Shower Opening Angle", Binning::Simple(44, -1.1, 1.1), kLeadingShowerCosAngle),
                              HistAxisVM("Leading Shower Momentum", Binning::Simple(100, 0, 5), kLeadingShowerEleP),
                              HistAxisVM("Truth Vertex Contained", Binning::Simple(2, 0, 2), kTFiducial),
                              HistAxisVM("Reco Shower Length", Binning::Simple(40, 0, 400), kRShowerLength),
                              HistAxisVM("Reco Shower End X", Binning::Simple(80, -400, 400), kRShowerEndX),
                              HistAxisVM("Reco Shower End Y", Binning::Simple(40, -200, 200), kRShowerEndY),
                              HistAxisVM("Reco Shower End Z", Binning::Simple(200, -1000, 1000), kRShowerEndZ),
                              HistAxisVM("Reco Shower Start X", Binning::Simple(80, -400, 400), kRShowerStartX),
                              HistAxisVM("Reco Shower Start Y", Binning::Simple(40, -200, 200), kRShowerStartY),
                              HistAxisVM("Reco Shower Start Z", Binning::Simple(200, -1000, 1000), kRShowerStartZ),
                              HistAxisVM("Cryostat Containing Vertex", Binning::Simple(2, 0, 2), kVertexCryo),
                              HistAxisVM("Truth PDG ID", Binning::Simple(6000, -3000, 3000), kTruthPDGID),
                              HistAxisVM("Is Clear Cosmic", Binning::Simple(2, 0, 2), rClearCosmic),
                              HistAxisVM("Largest Shower Passes Cut", Binning::Simple(2, 0, 2), rLargestShowerCut),
                              HistAxisVM("Good #{nu}_{e}", Binning::Simple(2, 0, 2), tGoodNue),
                              //HistAxisVM("All Shower Energies (GeV)", Binning::Simple(60,0.f,6.f),  kShowerE)
                              };

  std::vector<std::string> truthVarLabels = {"Number of True #nu_e", "Number of True #nu_e", "Number of True #nu_{#mu}"};
  std::vector<Binning> truthBins= {Binning::Simple(10, -5, 5), Binning::Simple(10, -5, 5), Binning::Simple(10, -5, 5)};
  std::vector<TruthVar> truthVars = {kTGoodNue, kTnNue, kTnNumu};
  std::vector<std::string> truthVarNames = {"NumNuE", "nNue", "nNumu"};

  std::vector<Var> treeVars = {kTruth_NeutrinoE, kLeadingShowerE, kLeadingShowerDedx, kLeadingShowerDensity, kLeadingShowerConversionGap,
                              kLongestTrack, kBaryDeltaZT, kBaryRadiusT, kRFiducial, kShwContained, kMuonChi2, kInteractionMode, kNeutrinoType, kSubLeadingShowerE,
                              kShowerEResidual, kSubleadingShowerEResidual, kNuScore, kTShowerLength, kShowerLenResidual, kVertexResidual,
                              kTShowerE, kTSubleadingShowerE, kTrackLenResidual, kTTrackLen, kLeadingShowerCosAngle, kLeadingShowerEleP,
                              kTFiducial, kRShowerLength, kRShowerEndX, kRShowerEndY, kRShowerEndZ, kRShowerStartX, kRShowerStartY, kRShowerStartZ,
                              kVertexCryo, kTruthPDGID, rClearCosmic, rLargestShowerCut, tGoodNue
                              };

  if ( useCVN ) {
    treeVars.insert(treeVars.end(), {kCVNNueScore, kCVNNumuScore, kCVNNcScore, kCVNCosmicScore});
    ax.insert(ax.end(),                              
        {HistAxisVM("CVN Nue Score", Binning::Simple(50, 0, 1), kCVNNueScore),
        HistAxisVM("CVN Numu Score", Binning::Simple(50, 0, 1), kCVNNumuScore),
        HistAxisVM("CVN NC Score", Binning::Simple(50, 0, 1), kCVNNcScore),
        HistAxisVM("CVN Cosmic Score", Binning::Simple(50, 0, 1), kCVNCosmicScore)}
    );
    varNames.insert(varNames.end(), {"kCVNNueScore", "kCVNNumuScore", "kCVNNcScore", "kCVNCosmicScore" });
  }
  
  std::vector<MultiVar> treeMultiVars;
  for (auto &v : treeVars) {
    treeMultiVars.emplace_back(
      ana::_MultiVar<caf::Proxy<caf::SRSlice>>{
        [v](const caf::Proxy<caf::SRSlice>* s) -> std::vector<double> {
          // wrap single Var result into a 1-element vector
          return { v(s) };
        }
      }
    );
  }  

  //treeMultiVars.emplace_back(kShowerE);
  //std::vector<std::string> treeVarNames = {"rAllShowerEnergies"};

  Cut treeCut = initialCut;
  for (auto& c : eShowerCuts) treeCut = treeCut && c;

  //SpillCut defaultSpillCut = kCRTPMT_correction;
  SpillCut defaultSpillCut = kNoSpillCut;

  const Tree t_sel ("t_sel", varNames, loader, treeMultiVars, defaultSpillCut, treeCut, kNoShift, true, true);
  const Tree t_eff ("t_eff", varNames, loader, treeMultiVars, defaultSpillCut, kSliceNeutrinoMatched, kNoShift, true, true);
  const Tree t_selEff ("t_selEff", varNames, loader, treeMultiVars, defaultSpillCut , kSliceNeutrinoMatched && treeCut, kNoShift, true, true);
  const Tree t_raw ("t_raw", varNames, loader, treeMultiVars, defaultSpillCut, kNoCut, kNoShift, true, true);
  const Tree t_evt ("t_evt", truthVarNames, loader, truthVars, defaultSpillCut, kNoTruthCut, kNoCut, kNoShift, true, true);

  // Add 2d histogram names after filling trees
  varNames.push_back("vertexCM");
  ax.push_back(HistAxisVM("Truth Vertex Contained", Binning::Simple(2, 0, 2), kTFiducial, "Reco Vertex Contained", Binning::Simple(2, 0, 2), kRFiducial));

  //eShower.createNM1Spectra(loader, axtnue);
  all.createNM1Spectra(loader, ax, varNames, defaultSpillCut);
  all.createSeqSpectra(loader, ax, varNames, defaultSpillCut);

  nue.createNM1Spectra(loader, ax, varNames, defaultSpillCut);
  nue.createSeqSpectra(loader, ax, varNames, defaultSpillCut);

  numu.createNM1Spectra(loader, ax, varNames, defaultSpillCut);
  numu.createSeqSpectra(loader, ax, varNames, defaultSpillCut);

  nc.createNM1Spectra(loader, ax, varNames, defaultSpillCut);
  nc.createSeqSpectra(loader, ax, varNames, defaultSpillCut);

  cosmic.createNM1Spectra(loader, ax, varNames, defaultSpillCut);
  cosmic.createSeqSpectra(loader, ax, varNames, defaultSpillCut);

  // kCacheMaxShowerSlice finds the slice ID with the max shower energy, only that slice is considered for the tNue spectra
  nueEff.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueEff.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);

  nueP.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueP.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);  

  nueA.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueA.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);

  fid.createNM1Spectra(loader, truthVarLabels, truthBins, truthVars, truthVarNames, kTGoodNueCut, kNoSpillCut);
  fid.createSeqSpectra(loader, truthVarLabels, truthBins, truthVars, truthVarNames, kTGoodNueCut, kNoSpillCut);

  nueCCQE.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueCCQE.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);

  nueCCRes.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueCCRes.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);

  nueCCDis.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueCCDis.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  
  nueCCCoh.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueCCCoh.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  
  nueCCMEC.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueCCMEC.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  
  nueNC.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  nueNC.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue);
  
  nueCCOther.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nueCCOther.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );

  nueEle.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nueEle.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );

  nueMu.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nueMu.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );

  nuePro.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nuePro.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );

  nuePi0.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nuePi0.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );

  nuePiC.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nuePiC.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );

  nuePh.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nuePh.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );

  nueOther.createNM1Spectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  nueOther.createSeqSpectra(loader, ax, varNames, defaultSpillCut && kTruth_goodNue );
  
  /////////////////////////
  // actually make the spectra
  loader.Go();  
  /////////////////////////////////////
 
  all.createNM1Hist();
  all.createSeqHist();

  nue.createNM1Hist();
  nue.createSeqHist();

  numu.createNM1Hist();
  numu.createSeqHist();

  nc.createNM1Hist();
  nc.createSeqHist();

  cosmic.createNM1Hist();
  cosmic.createSeqHist();

  fid.createNM1Hist();
  fid.createSeqHist();

  nueEff.createNM1Hist();
  nueEff.createSeqHist();

  nueP.createNM1Hist();
  nueP.createSeqHist();
  nueA.createNM1Hist();
  nueA.createSeqHist();

  nueCCQE.createNM1Hist();
  nueCCQE.createSeqHist();
  
  nueCCRes.createNM1Hist();
  nueCCRes.createSeqHist();
  
  nueCCDis.createNM1Hist();
  nueCCDis.createSeqHist();
  
  nueCCCoh.createNM1Hist();
  nueCCCoh.createSeqHist();
  
  nueCCMEC.createNM1Hist();
  nueCCMEC.createSeqHist();
  
  nueNC.createNM1Hist();
  nueNC.createSeqHist();
  
  nueCCOther.createNM1Hist();
  nueCCOther.createSeqHist();

  nueEle.createNM1Hist();
  nueEle.createSeqHist();

  nueMu.createNM1Hist();
  nueMu.createSeqHist();

  nuePro.createNM1Hist();
  nuePro.createSeqHist();

  nuePi0.createNM1Hist();
  nuePi0.createSeqHist();

  nuePiC.createNM1Hist();
  nuePiC.createSeqHist();

  nuePh.createNM1Hist();
  nuePh.createSeqHist();

  nueOther.createNM1Hist();
  nueOther.createSeqHist();

  ///////////////////
  // Save to file and make pretty plots
  ///////////////////
  f.cd();

  gDirectory->mkdir("allEvents");
  gDirectory->cd("allEvents");
  all.saveNM1Plots();
  all.saveSeqPlots();
  all.createCombinedNM1();
  all.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nue");
  gDirectory->cd("nue");
  nue.saveNM1Plots();
  nue.saveSeqPlots();
  nue.createCombinedNM1();
  nue.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("numu");
  gDirectory->cd("numu");
  numu.saveNM1Plots();
  numu.saveSeqPlots();
  numu.createCombinedNM1();
  numu.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nc");
  gDirectory->cd("nc");
  nc.saveNM1Plots();
  nc.saveSeqPlots();
  nc.createCombinedNM1();
  nc.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("cosmic");
  gDirectory->cd("cosmic");
  cosmic.saveNM1Plots();
  cosmic.saveSeqPlots();
  cosmic.createCombinedNM1();
  cosmic.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("fid");
  gDirectory->cd("fid");
  fid.saveNM1Plots();
  fid.saveSeqPlots();
  fid.createCombinedNM1();
  fid.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueEff");
  gDirectory->cd("nueEff");
  nueEff.saveNM1Plots();
  nueEff.saveSeqPlots();
  nueEff.createCombinedNM1();
  nueEff.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueP");;
  gDirectory->cd("nueP");
  nueP.saveNM1Plots();
  nueP.saveSeqPlots();
  nueP.createCombinedNM1();
  nueP.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueA");;
  gDirectory->cd("nueA");
  nueA.saveNM1Plots();  
  nueA.saveSeqPlots();
  nueA.createCombinedNM1();
  nueA.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueCCQE");
  gDirectory->cd("nueCCQE");
  nueCCQE.saveNM1Plots();
  nueCCQE.saveSeqPlots();
  nueCCQE.createCombinedNM1();
  nueCCQE.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueCCRes");
  gDirectory->cd("nueCCRes");
  nueCCRes.saveNM1Plots();
  nueCCRes.saveSeqPlots();
  nueCCRes.createCombinedNM1();
  nueCCRes.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueCCDis");
  gDirectory->cd("nueCCDis");
  nueCCDis.saveNM1Plots();
  nueCCDis.saveSeqPlots();
  nueCCDis.createCombinedNM1();
  nueCCDis.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueCCCoh");
  gDirectory->cd("nueCCCoh");
  nueCCCoh.saveNM1Plots();
  nueCCCoh.saveSeqPlots();
  nueCCCoh.createCombinedNM1();
  nueCCCoh.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueCCMEC");
  gDirectory->cd("nueCCMEC");
  nueCCMEC.saveNM1Plots();
  nueCCMEC.saveSeqPlots();
  nueCCMEC.createCombinedNM1();
  nueCCMEC.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueNC");
  gDirectory->cd("nueNC");
  nueNC.saveNM1Plots(); 
  nueNC.saveSeqPlots();
  nueNC.createCombinedNM1();
  nueNC.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueCCOther");
  gDirectory->cd("nueCCOther");
  nueCCOther.saveNM1Plots();
  nueCCOther.saveSeqPlots();
  nueCCOther.createCombinedNM1();
  nueCCOther.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueEle");
  gDirectory->cd("nueEle");
  nueEle.saveNM1Plots();
  nueEle.saveSeqPlots();
  nueEle.createCombinedNM1();
  nueEle.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueMu");
  gDirectory->cd("nueMu");
  nueMu.saveNM1Plots();
  nueMu.saveSeqPlots();
  nueMu.createCombinedNM1();
  nueMu.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nuePro");
  gDirectory->cd("nuePro");
  nuePro.saveNM1Plots();
  nuePro.saveSeqPlots();
  nuePro.createCombinedNM1();
  nuePro.createCombinedSeq();
  gDirectory->cd("../");  

  gDirectory->mkdir("nuePi0");
  gDirectory->cd("nuePi0");
  nuePi0.saveNM1Plots();
  nuePi0.saveSeqPlots();
  nuePi0.createCombinedNM1();
  nuePi0.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nuePiC");
  gDirectory->cd("nuePiC");
  nuePiC.saveNM1Plots();
  nuePiC.saveSeqPlots();
  nuePiC.createCombinedNM1();
  nuePiC.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nuePh");;
  gDirectory->cd("nuePh");
  nuePh.saveNM1Plots();
  nuePh.saveSeqPlots();
  nuePh.createCombinedNM1();
  nuePh.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("nueOther");
  gDirectory->cd("nueOther");
  nueOther.saveNM1Plots();
  nueOther.saveSeqPlots();
  nueOther.createCombinedNM1();
  nueOther.createCombinedSeq();
  gDirectory->cd("../");

  t_sel.SaveTo(gDirectory->mkdir("t_sel"));
  t_eff.SaveTo(gDirectory->mkdir("t_eff"));
  t_raw.SaveTo(gDirectory->mkdir("t_raw"));
  t_selEff.SaveTo(gDirectory->mkdir("t_selEff"));
  t_evt.SaveTo(gDirectory->mkdir("t_evt"));

  // gDirectory->mkdir("recoCutflows");
  // gDirectory->cd("recoCutflows");
  // std::vector<cutFlow<Cut>*> flows = {&nue, &numu, &nc, &cosmic};
  // std::vector<std::string> labels = {"#nu_e", "#nu_#mu", "NC", "Cosmic"};
  // compareCutflows<Cut>(f, flows, labels, "recoCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDedx", "rShowerDeDx", "recoCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerEnergy", "rShowerEnergy", "recoCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDensity", "rShowerDensity", "recoCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerConversionGap", "rShowerConversionGap", "recoCutflows");
  // gDirectory->cd("../");


  // gDirectory->mkdir("truthCutflows");
  // gDirectory->cd("truthCutflows");
  // flows = {&nue, &nueEff};
  // labels = {"#nu_e", "Truth #nu_e"};
  // compareCutflows<Cut>(f, flows, labels, "truthCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDedx", "rShowerDeDx", "truthCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerEnergy", "rShowerEnergy", "truthCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDensity", "rShowerDensity", "truthCutflows");
  // compareCutflowsNM1<Cut>(f, flows, labels, "ShowerConversionGap", "rShowerConversionGap", "truthCutflows");
  // gDirectory->cd("../");

  f.Close();
  
}  
