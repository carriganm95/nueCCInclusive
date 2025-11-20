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

void pionSelection(bool debug=false, std::string jobId="", std::string outputFile="debug.root"){

  std::cout << "jobId: " << jobId << " , " << std::endl;
  // /pnfs/sbn/data/sbn_fd/poms_production/data/Run2reprocess/reconstructed/icaruscode_v09_89_01_02p02/numimajority/flatcaf_prescaled/ //data
  TString fphase2 = "/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/00/*/*.flat.caf*.root";
  if ( debug ) fphase2 = "/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/05/04/detsim_2d_icarus_refactored_detsim_stage0_stage1_79315342_0.flat.caf-0985b172-de17-487e-8183-915eed082119.root";
  else if (jobId != "") fphase2 = jobId;

  std::cout << "Using input: " << fphase2 << std::endl;
  ///////////////////////
  //const std::string fphase2 = "/pnfs/sbn/data_add/sbn_nd/poms_production/official/MCP2025A/v10_04_03/prodoverlay_corsika_cosmics_proton_genie_rockbox_sce/LArv10/caf/05/0b/caf.flat.caf-3131a89e-10e8-4f4b-b738-af61313752e9.root";
  
  // Source of events
  SpectrumLoader loader(fphase2.Data());
  
  double POT = 3.0E20; //BNB 6.6 NuMI 6
  const double nominalIntensity =  5.13E13; // POT per spill, 5.13 one data (run 8515) when we get a stable beam 5E13 //tested with 6E13
  const double nominalLivetime = POT/nominalIntensity;

  TFile f(TString(outputFile), "RECREATE");

  const Cut totalCut = kNoCut;

  //std::vector<Cut> pionCuts = {totalCut, kRFiducialCut, kBaryDeltaZTCut, kMuonTrk, kMuonDeltaTrkVertexCut, kMuonCollectionHitsCut, kMuonTrkContainedCut, kMuonChi2Cut, kMuonTrkLenCut, goodPhotonCandidates, kShowersContained, primaryPhotonEnergyCut,  secondaryPhotonEnergyCut, nShowersCut};
  //std::vector<string> cutNames = {"NoCut", "Fiducial", "BarycenterDeltaZ", "MuonTrack", "MuonDeltaTrkVertex", "MuonCollectionHits", "MuonTrkContained", "MuonChi2", "MuonTrkLen", "PhotonCandidates", "ShowersContained", "primaryPhotonEnergy", "secondaryPhotonEnergy", "nShowers"};
  
  std::vector<Cut> pionCuts = {totalCut, kRFiducialCut, chargedPionVeto, kBarycenterCuts, kMuonTrk, kMuonDeltaTrkVertexCut, kMuonCollectionHitsCut, kMuonChi2Cut, kMuonTrkLenCut, goodPhotonCandidates, kShowersContained, primaryPhotonEnergyCut,  secondaryPhotonEnergyCut, nShowersCut};
  std::vector<string> cutNames = {"NoCut", "Fiducial", "chargedPionVeto", "Barycenter", "MuonTrack", "MuonDeltaTrkVertex", "MuonCollectionHits", "MuonChi2", "MuonTrkLen", "PhotonCandidates", "ShowersContained", "primaryPhotonEnergy", "secondaryPhotonEnergy", "nShowers"};
  cutFlow<Cut> all("all", pionCuts, cutNames);
  cutFlow<Cut> s1Mu1Pi0("s1Mu1Pi0", pionCuts, cutNames);
  cutFlow<Cut> b1MuNPi0("b1MuNPi0", pionCuts, cutNames);
  cutFlow<Cut> bNuMu("bNuMu", pionCuts, cutNames);
  cutFlow<Cut> bNuE("bNuE", pionCuts, cutNames);
  cutFlow<Cut> bNC("bNC", pionCuts, cutNames);
  cutFlow<Cut> bCosmic("bCosmic", pionCuts, cutNames);
  cutFlow<Cut> fid("fid", pionCuts, cutNames);


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

  std::vector<HistAxis> ax = {HistAxis("True E_{#nu} (GeV)", Binning::Simple(60,0.f,6.f),  kTruth_NeutrinoE),
                              HistAxis("Reco Leading Shower E (GeV)", Binning::Simple(20,0.f,1.f),  kLeadingShowerE),
                              HistAxis("Reco Subleading Shower E (GeV)", Binning::Simple(20, 0.f, 1.f), kSubLeadingShowerE),
                              HistAxis("Reco Leading Shower DeDx (MeV/cm)", Binning::Simple(50, 0.f, 10.f), kLeadingShowerDedx),
                              HistAxis("Reco Leading Shower Density (MeV/cm)", Binning::Simple(50, 0.f, 20.f), kLeadingShowerDensity),
                              HistAxis("Reco Leading Shower Conversion Gap (cm)", Binning::Simple(50, 0.f, 10.f), kLeadingShowerConversionGap),
                              HistAxis("Longest Track Length (cm)", Binning::Simple(40, 0, 400), kLongestTrack),
                              HistAxis("Barycenter Delta(Z) (cm)", Binning::Simple(200, 0, 200), kBaryDeltaZT),
                              HistAxis("Barycenter Radius (cm)", Binning::Simple(200, 0, 200), kBaryRadiusT),
                              HistAxis("Vertex in Fiducial Volume", Binning::Simple(2, 0, 2), kRFiducial),
                              HistAxis("Shower Contained in Fiducial Volume", Binning::Simple(2, 0, 2), kShwContained),
                              HistAxis("Muon Chi2 Score", Binning::Simple(50, 0, 200), kMuonChi2),
                              HistAxis("Cos(#theta(leading, sub-leading))", Binning::Simple(20, -1, 1), kShowerAngle),
                              HistAxis("#pi^{0} Mass (MeV)", Binning::Simple(20, 0, 200), kShowerInvariantMass),
                              HistAxis("Leading Shower Energy (GeV)", 50, 0, 5, kLeadingShowerE,
                                      "Subleading Shower Energy (GeV)", 50, 0, 5, kSubLeadingShowerE),
                              HistAxis("Number of Showers", 10, 0, 10, kNShowers),
                              HistAxis("Proton Chi2 Score", Binning::Simple(100, 0, 400), kProtonChi2),
                              HistAxis("#Delta(trk, vtx)", Binning::Simple(50, 0, 50), kMuonDeltaTrkVertex),
                              HistAxis("Collection Plane Hits", Binning::Simple(100, 0, 100), kMuonCollectionHits),
                              HistAxis("Trk Contained", Binning::Simple(3, -1, 2), kMuonTrkContained),
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
                    kCRTPMT_correction, initialCut && kNumuCC && kT1Mu1Pi0, kNoShift, true, true);

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

  t1Mu1Pi0.SaveTo(gDirectory->mkdir("signal"));
  allEvents.SaveTo(gDirectory->mkdir("all"));

  /*gDirectory->mkdir("recoCutflows");
  gDirectory->cd("recoCutflows");
  std::vector<cutFlow<Cut>*> flows = {&nue, &numu, &nc, &cosmic};
  std::vector<std::string> labels = {"#nu_e", "#nu_#mu", "NC", "Cosmic"};
  compareCutflows<Cut>(f, flows, labels, "recoCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDedx", "rShowerDeDx", "recoCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerEnergy", "rShowerEnergy", "recoCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDensity", "rShowerDensity", "recoCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerConversionGap", "rShowerConversionGap", "recoCutflows");
  gDirectory->cd("../");


  gDirectory->mkdir("truthCutflows");
  gDirectory->cd("truthCutflows");
  flows = {&nue, &tNue};
  labels = {"#nu_e", "Truth #nu_e"};
  compareCutflows<Cut>(f, flows, labels, "truthCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDedx", "rShowerDeDx", "truthCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerEnergy", "rShowerEnergy", "truthCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerDensity", "rShowerDensity", "truthCutflows");
  compareCutflowsNM1<Cut>(f, flows, labels, "ShowerConversionGap", "rShowerConversionGap", "truthCutflows");
  gDirectory->cd("../");*/

  /*hGoodMuon->Write("goodMuons");

  hTrkLen_numu->Write("hTrkLen_numu");

  h_p_mom_goodLead->Write("h_p_mom_goodLead");
  h_p_mom_good->Write("h_p_mom_good");
  h_p_mom_lead->Write("h_p_mom_lead");
  h_p_mom->Write("h_p_mom");

  h_m_mom_goodLead->Write("h_m_mom_goodLead");
  h_m_mom_good->Write("h_m_mom_good");
  h_m_mom_lead->Write("h_m_mom_lead");
  h_m_mom->Write("h_m_mom");*/

}  
