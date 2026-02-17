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

void pionSelection(bool debug=false, std::string jobId="", std::string outputFile="debug3.root", int start=-1, int stop=-1){

  std::cout << "jobId: " << jobId << " , " << std::endl;
  // /pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/numimajority/flatcaf_prescaled/ //data
  TString fphase2 = "/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/00/*/*.flat.caf*.root";
  if ( debug ) fphase2 = "/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/05/04/detsim_2d_icarus_refactored_detsim_stage0_stage1_79315342_0.flat.caf-0985b172-de17-487e-8183-915eed082119.root";
  else if (jobId != "") fphase2 = jobId;

  std::cout << "Using input: " << fphase2 << std::endl;
  ///////////////////////
  //const std::string fphase2 = "/pnfs/sbn/data_add/sbn_nd/poms_production/official/MCP2025A/v10_04_03/prodoverlay_corsika_cosmics_proton_genie_rockbox_sce/LArv10/caf/05/0b/caf.flat.caf-3131a89e-10e8-4f4b-b738-af61313752e9.root";
  std::string inputStr = fphase2.Data();
  std::cout << "Input string: " << inputStr << std::endl;

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

  //std::vector<Cut> pionCuts = {totalCut, kRFiducialCut, kBaryDeltaZTCut, kMuonTrk, kMuonDeltaTrkVertexCut, kMuonCollectionHitsCut, kMuonTrkContainedCut, kMuonChi2Cut, kMuonTrkLenCut, goodPhotonCandidates, kShowersContained, primaryPhotonEnergyCut,  secondaryPhotonEnergyCut, nShowersCut};
  //std::vector<string> cutNames = {"NoCut", "Fiducial", "BarycenterDeltaZ", "MuonTrack", "MuonDeltaTrkVertex", "MuonCollectionHits", "MuonTrkContained", "MuonChi2", "MuonTrkLen", "PhotonCandidates", "ShowersContained", "primaryPhotonEnergy", "secondaryPhotonEnergy", "nShowers"};
  
  std::vector<Cut> pionCuts = {totalCut, kRFiducialCut, chargedPionVeto, kBarycenterCuts, kMuonTrk, kMuonDeltaTrkVertexCut, kMuonCollectionHitsCut, kMuonChi2Cut, kMuonTrkLenCut, goodPhotonCandidates, openAngle, kShowersContained, primaryPhotonEnergyCut,  secondaryPhotonEnergyCut};
  std::vector<string> cutNames = {"NoCut", "Fiducial", "chargedPionVeto", "Barycenter", "MuonTrack", "MuonDeltaTrkVertex", "MuonCollectionHits", "MuonChi2", "MuonTrkLen", "PhotonCandidates", "openAngle", "ShowersContained", "primaryPhotonEnergy", "secondaryPhotonEnergy"};
  cutFlow<Cut> all("all", pionCuts, cutNames);
  cutFlow<Cut> s1Mu1Pi0("s1Mu1Pi0", pionCuts, cutNames);
  cutFlow<Cut> b1MuNPi0("b1MuNPi0", pionCuts, cutNames);
  cutFlow<Cut> bNuMu("bNuMu", pionCuts, cutNames);
  cutFlow<Cut> bNuE("bNuE", pionCuts, cutNames);
  cutFlow<Cut> bNC("bNC", pionCuts, cutNames);
  cutFlow<Cut> bCosmic("bCosmic", pionCuts, cutNames);
  cutFlow<Cut> fid("fid", pionCuts, cutNames);
  
  cutFlow<Cut> allEle("allEle", pionCuts, cutNames); // all events (true electron only)
  cutFlow<Cut> allMu("allMu", pionCuts, cutNames); // all events (muon only)
  cutFlow<Cut> allPro("allPro", pionCuts, cutNames); // all events (proton only)
  cutFlow<Cut> allPh("allPh", pionCuts, cutNames); // all events (photon only)
  cutFlow<Cut> allPi0("allPi0", pionCuts, cutNames); // all events (neutral pion only)
  cutFlow<Cut> allPiC("allPiC", pionCuts, cutNames); // all events (charged pion only)
  cutFlow<Cut> allOther("allOther", pionCuts, cutNames); // all events (other only)

  Cut initialCut = kNotClearCosmic;

  all.createNM1Cuts(initialCut);
  all.createSeqCuts(initialCut);
  s1Mu1Pi0.createNM1Cuts(initialCut && kNumuCC && kT1Mu1Pi0);
  s1Mu1Pi0.createSeqCuts(initialCut && kNumuCC && kT1Mu1Pi0);
  b1MuNPi0.createNM1Cuts(initialCut && kNumuCC && kT1MuNPi0);
  b1MuNPi0.createSeqCuts(initialCut && kNumuCC && kT1MuNPi0);
  bNuMu.createNM1Cuts(initialCut && kNumuCC && !kT1Mu1Pi0 && !kT1MuNPi0);
  bNuMu.createSeqCuts(initialCut && kNumuCC && !kT1Mu1Pi0 && !kT1MuNPi0);
  bNuE.createNM1Cuts(initialCut && kNueCC);
  bNuE.createSeqCuts(initialCut && kNueCC);
  bNC.createNM1Cuts(initialCut && kNC);
  bNC.createSeqCuts(initialCut && kNC);
  bCosmic.createNM1Cuts(initialCut && kThisCosmic);
  bCosmic.createSeqCuts(initialCut && kThisCosmic);
  fid.createNM1Cuts(initialCut);
  fid.createSeqCuts(initialCut);
  
  allEle.createNM1Cuts(initialCut && kIsTrueElectron);
  allEle.createSeqCuts(initialCut && kIsTrueElectron);
  allMu.createNM1Cuts(initialCut && kIsTrueMuon);
  allMu.createSeqCuts(initialCut && kIsTrueMuon);
  allPro.createNM1Cuts(initialCut && kIsTrueProton);
  allPro.createSeqCuts(initialCut && kIsTrueProton);
  allPh.createNM1Cuts(initialCut && kIsTruePhoton);
  allPh.createSeqCuts(initialCut && kIsTruePhoton);
  allPi0.createNM1Cuts(initialCut && kIsTrueNeutralPion);
  allPi0.createSeqCuts(initialCut && kIsTrueNeutralPion);
  allPiC.createNM1Cuts(initialCut && kIsTrueChargedPion);
  allPiC.createSeqCuts(initialCut && kIsTrueChargedPion);
  allOther.createNM1Cuts(initialCut && kOtherTruthPDG);
  allOther.createSeqCuts(initialCut && kOtherTruthPDG);

  // Spectrum
  ////######################################################/////
 
  const HistAxis axtnue  ("True E_{#nu} (GeV)",    Binning::Simple(60,0.f,6.f),  kTruth_NeutrinoE); 
  const HistAxis ax_recoShowerE ("Reco Shower E (GeV)", Binning::Simple(250, 0, 5), kLeadingShowerE);
  const HistAxis ax_leadingProtonP ("Leading Proton Momentum (GeV)", Binning::Simple(200, 0, 2), kLeadingProtonMomentum);
  const HistAxis ax_leadingMuonP ("Leading Muon Momentum (GeV)", Binning::Simple(200, 0, 2), kLeadingMuonMomentum);
  const HistAxis ax_fiducialShower("Shower Contained in Fiducial Volume", Binning::Simple(2, 0, 2), kShwContained);

  std::vector<std::string> varNames = {"tNueEnergy", "rShowerEnergy", "rSubleadingShowerEnergy", "rShowerDeDx", "rShowerDensity", "rShowerConversionGap", "rLongestTrackLength", 
                                      "rBarycenterDeltaZ", "rBarycenterRadius", "rFiducial", "rShowerContained", "rMuonChi2", "PhotonAngle", "Pi0Mass", "ShowerEnergies", "numShowers",
                                      "rProtonChi2", "rMuonDeltaTrkVertex", "rMuonCollectionHits", "rMuonTrkContained"};

  std::vector<HistAxisVM> ax = {HistAxisVM("True E_{#nu} (GeV)", Binning::Simple(60,0.f,6.f),  kTruth_NeutrinoE),
                              HistAxisVM("Reco Leading Shower E (GeV)", Binning::Simple(20,0.f,1.f),  kLeadingShowerE),
                              HistAxisVM("Reco Subleading Shower E (GeV)", Binning::Simple(20, 0.f, 1.f), kSubLeadingShowerE),
                              HistAxisVM("Reco Leading Shower DeDx (MeV/cm)", Binning::Simple(50, 0.f, 10.f), kLeadingShowerDedx),
                              HistAxisVM("Reco Leading Shower Density (MeV/cm)", Binning::Simple(50, 0.f, 20.f), kLeadingShowerDensity),
                              HistAxisVM("Reco Leading Shower Conversion Gap (cm)", Binning::Simple(50, 0.f, 10.f), kLeadingShowerConversionGap),
                              HistAxisVM("Longest Track Length (cm)", Binning::Simple(40, 0, 400), kLongestTrack),
                              HistAxisVM("Barycenter Delta(Z) (cm)", Binning::Simple(200, 0, 200), kBaryDeltaZT),
                              HistAxisVM("Barycenter Radius (cm)", Binning::Simple(200, 0, 200), kBaryRadiusT),
                              HistAxisVM("Vertex in Fiducial Volume", Binning::Simple(2, 0, 2), kRFiducial),
                              HistAxisVM("Shower Contained in Fiducial Volume", Binning::Simple(2, 0, 2), kShwContained),
                              HistAxisVM("Muon Chi2 Score", Binning::Simple(50, 0, 200), kMuonChi2),
                              HistAxisVM("Cos(#theta(leading, sub-leading))", Binning::Simple(20, -1, 1), kShowerAngle),
                              HistAxisVM("#pi^{0} Mass (MeV)", Binning::Simple(50, 0, 500), kShowerInvariantMass),
                              HistAxisVM("Leading Shower Energy (GeV)", Binning::Simple(50, 0, 5), kLeadingShowerE,
                                      "Subleading Shower Energy (GeV)", Binning::Simple(50, 0, 5), kSubLeadingShowerE),
                              HistAxisVM("Number of Showers", Binning::Simple(10, 0, 10), kNShowers),
                              HistAxisVM("Proton Chi2 Score", Binning::Simple(100, 0, 400), kProtonChi2),
                              HistAxisVM("#Delta(trk, vtx)", Binning::Simple(50, 0, 50), kMuonDeltaTrkVertex),
                              HistAxisVM("Collection Plane Hits", Binning::Simple(100, 0, 100), kMuonCollectionHits),
                              HistAxisVM("Trk Contained", Binning::Simple(3, -1, 2), kMuonTrkContained),
                              };

  std::vector<std::string> truthVarLabels = {"Number of True 1Mu1Pi0"};
  std::vector<Binning> truthBins= {Binning::Simple(10, -5, 5)};
  std::vector<TruthVar> truthVars = {kTGood1Mu1Pi0};
  std::vector<std::string> truthVarNames = {"NumSignal"};


  const Tree t1Mu1Pi0 ("t1Mu1Pi0", 
                    {"tNeutrinoE", "rLeadingShowerE", "rSubLeadingShowerE", "rLeadingShowerDeDx",
                    "rLeadingShowerDensity", "rLeadingShowerConversionGap", "rLongestTrack", "rBarycenterDeltaZT",
                    "rBarycenterRadiusT", "rFiducial", "rShowerFiducial", "rMuonChi2", "rShowerCosTheta", "rShowerM", "rNShowers",
                    "rProtonChi2", "rDeltaTrkVtx", "rNCollectionPlaneHits", "rTrkContained"}, 
                    loader, 
                    {kTruth_NeutrinoE, kLeadingShowerE, kSubLeadingShowerE, kLeadingShowerDedx, 
                    kLeadingShowerDensity, kLeadingShowerConversionGap, kLongestTrack, kBaryDeltaZT, 
                    kBaryRadiusT, kRFiducial, kShwContained, kMuonChi2, kShowerAngle, kShowerInvariantMass,
                    kNShowers, kProtonChi2, kMuonDeltaTrkVertex, kMuonCollectionHits, kMuonTrkContained}, 
                    kCRTPMT_correction, initialCut && kNumuCC && kT1Mu1Pi0 && nShowersCut, kNoShift, true, true);

    const Tree allEvents ("allEvents", 
                    {"tNeutrinoE", "rLeadingShowerE", "rSubLeadingShowerE", "rLeadingShowerDeDx",
                    "rLeadingShowerDensity", "rLeadingShowerConversionGap", "rLongestTrack", "rBarycenterDeltaZT",
                    "rBarycenterRadiusT", "rFiducial", "rShowerFiducial", "rMuonChi2", "rShowerCosTheta", "rShowerM", "rNShowers",
                    "rProtonChi2", "rDeltaTrkVtx", "rNCollectionPlaneHits", "rTrkContained"}, 
                    loader, 
                    {kTruth_NeutrinoE, kLeadingShowerE, kSubLeadingShowerE, kLeadingShowerDedx, 
                    kLeadingShowerDensity, kLeadingShowerConversionGap, kLongestTrack, kBaryDeltaZT, 
                    kBaryRadiusT, kRFiducial, kShwContained, kMuonChi2, kShowerAngle, kShowerInvariantMass,
                    kNShowers, kProtonChi2, kMuonDeltaTrkVertex, kMuonCollectionHits, kMuonTrkContained}, 
                    kCRTPMT_correction, initialCut, kNoShift, true, true);

    const Tree signalShowerE0 ("signalShowerE0", 
      {"rShowerE0", 
        "tShowerE_E0",
        "tShowerE_W0"}, 
      loader, 
      {kShowerE0, 
        tShowerE_E0, 
        tShowerE_W0},  
      kCRTPMT_correction, 
      initialCut && kNumuCC && kT1Mu1Pi0 && nShowersCut, 
      kNoShift, 
      true, 
      true);
    const Tree signalShowerE1 ("signalShowerE1", 
      {"rShowerE1",
        "tShowerE_E1", 
        "tShowerE_W1"}, 
      loader, 
      {kShowerE1,
        tShowerE_E1, 
        tShowerE_W1},  
      kCRTPMT_correction, 
      initialCut && kNumuCC && kT1Mu1Pi0 && nShowersCut, 
      kNoShift, 
      true, 
      true);
    const Tree signalShowerE2 ("signalShowerE2", 
      {"rShowerE2", 
        "tShowerE_E2",
        "tShowerE_W2"}, 
      loader, 
      {kShowerE2, 
        tShowerE_E2, 
        tShowerE_W2},  
      kCRTPMT_correction, 
      initialCut && kNumuCC && kT1Mu1Pi0 && nShowersCut, 
      kNoShift, 
      true, 
      true);

  all.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  all.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction);  
  s1Mu1Pi0.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  s1Mu1Pi0.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  b1MuNPi0.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  b1MuNPi0.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  bNuMu.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  bNuMu.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  bNuE.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  bNuE.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  bNC.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  bNC.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  bCosmic.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  bCosmic.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  fid.createNM1Spectra(loader, truthVarLabels, truthBins, truthVars, truthVarNames, kTGood1Mu1Pi0Cut, kNoSpillCut);
  fid.createSeqSpectra(loader, truthVarLabels, truthBins, truthVars, truthVarNames, kTGood1Mu1Pi0Cut, kNoSpillCut);

  allEle.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  allEle.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction);  
  allMu.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  allMu.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  allPro.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  allPro.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  allPh.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  allPh.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  allPi0.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  allPi0.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  allPiC.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  allPiC.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction); 
  allOther.createNM1Spectra(loader, ax, varNames, kCRTPMT_correction);
  allOther.createSeqSpectra(loader, ax, varNames, kCRTPMT_correction);

  // neutrino energy truth source
  Spectrum s_energy_1mu1pi0 (loader, ax_fiducialShower, kCRTPMT_correction, initialCut && kNumuCC && kT1Mu1Pi0);
  Spectrum s_energy_Nmu1pi0 (loader, ax_fiducialShower, kCRTPMT_correction, initialCut && kNumuCC && kT1MuNPi0);
  Spectrum s_energy_numuOther (loader, ax_fiducialShower, kCRTPMT_correction, initialCut && kNumuCC && !kT1Mu1Pi0 && !kT1MuNPi0);
  Spectrum s_energy_nue (loader,  ax_fiducialShower, kCRTPMT_correction, initialCut && kNueCC); // energy of true nue
  Spectrum s_energy_other(loader,  ax_fiducialShower, kCRTPMT_correction, initialCut && kNC); // energy of true nc nu
  Spectrum s_energy_cos  (loader,  ax_fiducialShower, kCRTPMT_correction, initialCut && kThisCosmic); // energy of true cosmic nu
  Spectrum s_energy_all (loader, ax_fiducialShower, kCRTPMT_correction, initialCut);
  
  /////////////////////////
  // actually make the spectra
  loader.Go();  
  /////////////////////////////////////
 
  TH1* h_energy_1mu1pi0 = s_energy_1mu1pi0.ToTH1(POT);
  TH1* h_energy_Nmu1pi0 = s_energy_Nmu1pi0.ToTH1(POT);
  TH1* h_energy_numu = s_energy_numuOther.ToTH1(POT);
  TH1* h_energy_nue = s_energy_nue.ToTH1(POT);
  TH1* h_energy_other = s_energy_other.ToTH1(POT);
  TH1* h_energy_cos = s_energy_cos.ToTH1(POT);
  TH1* h_energy_all = s_energy_all.ToTH1(POT);

  all.createNM1Hist();
  all.createSeqHist();
  s1Mu1Pi0.createNM1Hist();
  s1Mu1Pi0.createSeqHist();
  b1MuNPi0.createNM1Hist();
  b1MuNPi0.createSeqHist();
  bNuMu.createNM1Hist();
  bNuMu.createSeqHist();
  bNuE.createNM1Hist();
  bNuE.createSeqHist();
  bNC.createNM1Hist();
  bNC.createSeqHist();
  bCosmic.createNM1Hist();
  bCosmic.createSeqHist();
  fid.createNM1Hist();
  fid.createSeqHist();

  allEle.createNM1Hist();
  allEle.createSeqHist();
  allMu.createNM1Hist();
  allMu.createSeqHist();
  allPro.createNM1Hist();
  allPro.createSeqHist();
  allPh.createNM1Hist();
  allPh.createSeqHist();
  allPi0.createNM1Hist();
  allPi0.createSeqHist();
  allPiC.createNM1Hist();
  allPiC.createSeqHist();
  allOther.createNM1Hist();
  allOther.createSeqHist();

  ///////////////////
  // Save to file and make pretty plots
  ///////////////////
  f.cd();
  h_energy_1mu1pi0->Write("h_energy_1mu1pi0");
  h_energy_Nmu1pi0->Write("h_energy_Nmu1pi0");
  h_energy_numu -> Write("h_energy_numu");
  h_energy_other -> Write("h_energy_nc");
  h_energy_nue -> Write("h_energy_nue");
  h_energy_cos -> Write("h_energy_cosmic");
  h_energy_all->Write("h_energy_all");

  gDirectory->mkdir("allEvents");
  gDirectory->cd("allEvents");
  all.saveNM1Plots();
  all.saveSeqPlots();
  all.createCombinedNM1();
  all.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("s1Mu1Pi0");
  gDirectory->cd("s1Mu1Pi0");
  s1Mu1Pi0.saveNM1Plots();
  s1Mu1Pi0.saveSeqPlots();
  s1Mu1Pi0.createCombinedNM1();
  s1Mu1Pi0.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("b1MuNPi0");
  gDirectory->cd("b1MuNPi0");
  b1MuNPi0.saveNM1Plots();
  b1MuNPi0.saveSeqPlots();
  b1MuNPi0.createCombinedNM1();
  b1MuNPi0.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("bNuMu");
  gDirectory->cd("bNuMu");
  bNuMu.saveNM1Plots();
  bNuMu.saveSeqPlots();
  bNuMu.createCombinedNM1();
  bNuMu.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("bNuE");
  gDirectory->cd("bNuE");
  bNuE.saveNM1Plots();
  bNuE.saveSeqPlots();
  bNuE.createCombinedNM1();
  bNuE.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("bNC");
  gDirectory->cd("bNC");
  bNC.saveNM1Plots();
  bNC.saveSeqPlots();
  bNC.createCombinedNM1();
  bNC.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("bCosmic");
  gDirectory->cd("bCosmic");
  bCosmic.saveNM1Plots();
  bCosmic.saveSeqPlots();
  bCosmic.createCombinedNM1();
  bCosmic.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("fid");
  gDirectory->cd("fid");
  fid.saveNM1Plots();
  fid.saveSeqPlots();
  fid.createCombinedNM1();
  fid.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("allEle");
  gDirectory->cd("allEle");
  allEle.saveNM1Plots();
  allEle.saveSeqPlots();
  allEle.createCombinedNM1();
  allEle.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("allMu");
  gDirectory->cd("allMu");
  allMu.saveNM1Plots();
  allMu.saveSeqPlots();
  allMu.createCombinedNM1();
  allMu.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("allPro");;
  gDirectory->cd("allPro");
  allPro.saveNM1Plots();
  allPro.saveSeqPlots();
  allPro.createCombinedNM1();
  allPro.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("allPh");;
  gDirectory->cd("allPh");
  allPh.saveNM1Plots();
  allPh.saveSeqPlots();
  allPh.createCombinedNM1();
  allPh.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("allPi0");;
  gDirectory->cd("allPi0");
  allPi0.saveNM1Plots();
  allPi0.saveSeqPlots();
  allPi0.createCombinedNM1();
  allPi0.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("allPiC");;
  gDirectory->cd("allPiC");
  allPiC.saveNM1Plots();
  allPiC.saveSeqPlots();
  allPiC.createCombinedNM1();
  allPiC.createCombinedSeq();
  gDirectory->cd("../");

  gDirectory->mkdir("allOther");;
  gDirectory->cd("allOther");
  allOther.saveNM1Plots();
  allOther.saveSeqPlots();
  allOther.createCombinedNM1();
  allOther.createCombinedSeq();
  gDirectory->cd("../");

  t1Mu1Pi0.SaveTo(gDirectory->mkdir("signal"));
  allEvents.SaveTo(gDirectory->mkdir("all"));
  signalShowerE0.SaveTo(gDirectory->mkdir("signalShowerE0"));
  signalShowerE1.SaveTo(gDirectory->mkdir("signalShowerE1"));
  signalShowerE2.SaveTo(gDirectory->mkdir("signalShowerE2"));

}  
