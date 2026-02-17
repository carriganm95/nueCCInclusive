#include "nuCCInclusiveCuts.h"

// The below is a strange workaround to pass spill information to a slice cut
// This is needed to ensure that only one slice per spill is marked as containing the true neutrino
// If the functions are not used for the spillCut and Cut there is a library linking issue, doesn't cause a crash but is ugly (instead make the code ugly..)
// These functions are constructed at first use rather than at module load time -> avoid some linkage issues
namespace {
  // File-local singletons via function-local statics
  const caf::SRSpillProxy*& SpillCache()
  {
    static const caf::SRSpillProxy* spill = nullptr;
    return spill;
  }

  int& MaxIdxCache()
  {
    static int idx = -1;
    return idx;
  }

  int FindMaxShowerSliceIdx(const caf::SRSpillProxy* sr)
  {
    int bestIdx = -1;
    double maxE = -1.0;
    for (unsigned i = 0; i < sr->slc.size(); ++i) {
      for (auto& pfp : sr->slc[i].reco.pfp) {
        if (!ana::checkShower(pfp, &sr->slc[i])) continue;
        const double e = pfp.shw.plane[2].energy;
        if (e > maxE) { maxE = e; bestIdx = (int)i; }
      }
    }
    return bestIdx;
  }
}

// Important: return references to function-static Cuts
static const ana::SpillCut& kCacheMaxShowerSlice()
{
  static const ana::SpillCut cut([](const caf::SRSpillProxy* sr) -> bool {
    SpillCache()  = sr;
    MaxIdxCache() = FindMaxShowerSliceIdx(sr);
    return MaxIdxCache() >= 0;  // Or return true if you don't want to reject spills
  });
  return cut;
}

static const ana::Cut& kIsMaxShowerEnergySlice()
{
  static const ana::Cut cut([](const caf::SRSliceProxy* slc) -> bool {
    auto* sr = SpillCache();
    if (!sr) return false;
    const int idx = MaxIdxCache();
    if (idx < 0 || (size_t)idx >= sr->slc.size()) return false;
    return (&sr->slc[idx] == slc);
  });
  return cut;
}

//   // Helper function
// int FindMaxShowerSliceIdx(const caf::SRSpillProxy* sr) {
//   int bestIdx = -1;
//   double maxEnergy = -1;
  
//   for (unsigned int i = 0; i < sr->slc.size(); i++) {
//     for (auto& pfp : sr->slc[i].reco.pfp) {
//       if (!checkShower(pfp, &sr->slc[i])) continue;
//       if (pfp.shw.plane[2].energy > maxEnergy) {
//         maxEnergy = pfp.shw.plane[2].energy;
//         bestIdx = i;
//       }
//     }
//   }
//   return bestIdx;
// }

// //SpillCut that caches the spill and computes max slice index
// const SpillCut kCacheMaxShowerSlice([](const caf::SRSpillProxy* sr) -> bool {
//   g_CurrentSpill = sr;
//   g_MaxShowerSliceIdx = FindMaxShowerSliceIdx(sr);
//   return (g_MaxShowerSliceIdx >= 0);  // Reject spills with no valid slice
// });

// //Slice Cut that checks if this is the max shower slice
// const Cut kIsMaxShowerEnergySlice([](const caf::SRSliceProxy* slc) -> bool {
//   if (!g_CurrentSpill || g_MaxShowerSliceIdx < 0) return false;
//   return (&g_CurrentSpill->slc[g_MaxShowerSliceIdx] == slc);
// });

namespace ana {

static std::string trackScoreBranch = "trackScore";
void SetTrackScoreBranch(const std::string& branchName) {
  trackScoreBranch = branchName;
}

////////////////////////////////
///////// Helper Functions /////
////////////////////////////////

float trkDistance(caf::Proxy<caf::SRPFP>& p, const caf::SRSliceProxy* slc){
  float dist = -1;
  dist = sqrt( pow(p.trk.start.x - slc->vertex.x, 2) + pow(p.trk.start.y - slc->vertex.y, 2) + pow(p.trk.start.z - slc->vertex.z, 2) );
  //cout << "Distance  between track and vertex " << dist << std::endl;
  return dist;
}

bool trackFit(caf::Proxy<caf::SRPFP>& p){
  if ( isnan(p.trk.start.x) ) return false; // track start x cannot be Nan
  if ( ( p.trk.len <= 0 ) || ( isnan(p.trk.len) ) ) return false; //track length > 0 and not Nan
  if ( p.trk.calo[p.trk.bestplane].nhit < 5 ) return false; // > 5 hits in collection plane
  return true;
}

bool containedTrack(caf::Proxy<caf::SRPFP>& p, float dist=10){
  // if track endpoint within 10cm of boundaries check chi2
  if ( isnan(p.trk.end.x) || isnan(p.trk.end.y) || isnan(p.trk.end.z) ) return false; //make sure track endpoints are not Nan
  if ( !(p.trk.end.x < -61.94 - dist && p.trk.end.x > -358.49 + dist) && !(p.trk.end.x > 61.94 + dist && p.trk.end.x < 358.49 - dist) ) return false; // x endpoint contined
  if ( !(p.trk.end.y > -181.86 + dist && p.trk.end.y < 134.96 - dist) ) return false; // y endpoint contained
  if ( !(p.trk.end.z > -894.951 + dist && p.trk.end.z < 894.951 - dist) ) return false; // z endpoint contained
  return true;
}

bool containedVertex(const caf::Proxy<caf::SRVector3D>& v, float x_low=10, float x_high=10, float y_low=10, float y_high=10, float z_low=10, float z_high=10){
  // if track endpoint within 10cm of boundaries check chi2
  if ( isnan(v.x) || isnan(v.y) || isnan(v.z) ) return false; //make sure vertex points are not Nan
  if ( !(v.x < -61.94 - x_low && v.x > -358.49 + x_high) && !(v.x > 61.94 + x_high && v.x < 358.49 - x_low) ) return false; // x vertex contined
  if ( !(v.y > -181.86 + y_high && v.y < 134.96 - y_low) ) return false; // y vertex contained
  if ( !(v.z > -894.951 + z_high && v.z < 894.951 - z_low) ) return false; // z vertex contained
  return true;
}

bool containedShower(caf::Proxy<caf::SRPFP>& p){
  
  double this_endx = p.shw.start.x + (p.shw.dir.x * p.shw.len);
  double this_endy = p.shw.start.y + (p.shw.dir.y * p.shw.len);
  double this_endz = p.shw.start.z + (p.shw.dir.z * p.shw.len);

  bool startx = ((-358.49 + 5. < p.shw.start.x) && (p.shw.start.x < -61.94 - 5.))
                || ((61.94 + 5. < p.shw.start.x) && (p.shw.start.x < 358.49 - 5. ));
  bool endx   = ((-358.49 + 5. < this_endx) && (this_endx < -61.94 - 5.))
                || ((61.94 + 5. < this_endx) && (this_endx < 358.49 - 5.));

  bool starty = (-181.86 + 5. < p.shw.start.y) && (p.shw.start.y < 134.96 - 5.);
  bool endy   = (-181.86 + 5. < this_endy) && (this_endy < 134.96 - 5.);

  bool startz = (-894.951 + 5. < p.shw.start.z) && (p.shw.start.z < 894.951 - 5.);
  bool endz   = (-894.951 + 5. < this_endz) && (this_endz < 894.951 - 5.);

  return (startx && endx && starty && endy && startz && endz);
}

bool IsTracklikeTrack( caf::Proxy<caf::SRPFP>& p, int id = -1) {
  double score = -1;
  if (trackScoreBranch == "trackScore") {
    return (!std::isnan(p.trackScore) && p.trackScore > 0.45);
  }
  else if (trackScoreBranch == "ngscore") {
    return (!std::isnan(p.ngscore.sem_cat) && p.ngscore.sem_cat == id);
  }
  else {
    std::cerr << "Unknown track score branch: " << trackScoreBranch << std::endl;
    return false; // or throw an exception
  }
}

bool checkMuon(caf::Proxy<caf::SRPFP>& p, const caf::SRSliceProxy* slc){
  if ( !IsTracklikeTrack(p, 0) ) return false; //make sure track score is good
  if ( std::isnan(p.trk.start.x) || std::isnan(p.trk.len) || p.trk.len <= 0. ) return false; //track start and len must be valid
  if ( trkDistance(p, slc) >= 10 ) return false; // dist(track, vertex) < 10 cm
  if ( !p.parent_is_primary ) return false; // must be a primary particle
  if ( p.trk.calo[2].nhit <= 5 ) return false; //require at least 5 collection plane hits
  if ( containedTrack(p) && (p.trk.chi2pid[2].chi2_muon >= 30 || p.trk.chi2pid[2].chi2_proton <= 60) ) return false; // check if track is contained within 10cm, if so needs good muon/proton scores
  return true; 
}

bool checkProton(caf::Proxy<caf::SRPFP>& p, const caf::SRSliceProxy* slc){
  if ( !IsTracklikeTrack(p, 1) ) return false; //make sure track score is good
  if ( std::isnan(p.trk.start.x) || std::isnan(p.trk.len) || p.trk.len <= 0. ) return false; //track start and len must be valid
  if ( trkDistance(p, slc) > 10 ) return false; // dist(track, vertex) < 10 cm
  if ( !p.parent_is_primary ) return false; // must be a primary particle
  if ( p.trk.calo[2].nhit < 5 ) return false; //require at least 5 collection plane hits
  if ( containedTrack(p) && (p.trk.chi2pid[2].chi2_muon <= 30 || p.trk.chi2pid[2].chi2_proton >= 90) ) return false; // check if track is contained within 10cm, if so needs good muon/proton scores

  return true;
}

bool checkChargedPion(caf::Proxy<caf::SRPFP>& p, const caf::SRSliceProxy* slc){
  if ( std::isnan(p.trk.start.x) || std::isnan(p.trk.len) || p.trk.len <= 0. ) return false; //track start and len must be valid
  if ( trkDistance(p, slc) > 10 ) return false; // dist(track, vertex) < 10 cm
  if ( !p.parent_is_primary ) return false; // must be a primary particle
  if ( p.trk.calo[2].nhit < 5 ) return false; //require at least 5 collection plane hits
  if ( containedTrack(p) && (p.trk.chi2pid[2].chi2_muon >= 30 || p.trk.chi2pid[2].chi2_proton <= 60) ) return false; // check if track is contained within 10cm, if so needs good muon/proton scores

  return true;
}

bool checkShower(caf::Proxy<caf::SRPFP>& p, const caf::SRSliceProxy* slc, int plane){
  if ( trackScoreBranch == "trackScore"){
    if ( std::isnan(p.trackScore) || p.trackScore < 0 || p.trackScore > 0.45) return false;
  }
  else if ( trackScoreBranch == "ngscore"){
    if ( std::isnan(p.ngscore.sem_cat) || p.ngscore.sem_cat != 2 ) return false;
  }
  else {
    std::cerr << "Unknown track score branch: " << trackScoreBranch << std::endl;
    return false;
  }
  if ( std::isnan(p.shw.plane[plane].energy) || p.shw.plane[plane].energy < 0 ) return false;
  if ( !p.parent_is_primary ) return false;
  return true;
}

// Non-default overload to avoid ambiguity when calling with two arguments
bool checkShower(caf::Proxy<caf::SRPFP>& p, const caf::SRSliceProxy* slc){
  return checkShower(p, slc, 2);
}

bool checkEle(caf::Proxy<caf::SRPFP>& s, const caf::SRSliceProxy* slc, int Nm1 = -1){

  int cutIndex = 0;
  int longestTrackId = kLongestTrackId(slc);

  //cout << "vertex: " << slc->vertex.x << ", " << slc->vertex.y << ", " << slc->vertex.z << ", " << containedVertex(slc->vertex, 25, 25, 25, 25, 50, 30) << std::endl;
  //if ( !containedVertex(slc->vertex, 25, 25, 25, 25, 50, 30) ) return false; // vertex must be contained in detector volume
  //if ( longestTrackId >= 0 && slc->reco.pfp.at(longestTrackId).trk.len >= 110 ) return false; // longest track must be < 110 cm
  //if ( !containedShower(s) ) return false; // shower must be contained 5cm from edges

  if ( Nm1 != cutIndex && s.shw.plane[2].energy <= 0.2 ) return false; // shower energy > 200
  cutIndex++;
  if ( Nm1 != cutIndex && s.shw.plane[2].dEdx >= 3.625 ) return false; // largest shower dEdx < 3.625 MeV/cm
  cutIndex++;
  if ( Nm1 != cutIndex && s.shw.conversion_gap >= 3.25 ) return false; // conversion gap < 3.25 cm
  cutIndex++;
  if ( Nm1 != cutIndex && s.shw.density <= 4.5 ) return false; // shower density > 4.5 MeV/cm

  return true;
}

bool ShowerContained(const caf::SRSliceProxy* slc, int id=-1){

  if (id == -1) id = kLargestShowerId(slc);
  if ( id == -1 ) return false;

  ///slc->rceo.pfp[largestShwIdx].shw
  double this_endx = slc->reco.pfp[id].shw.start.x + (slc->reco.pfp[id].shw.dir.x * slc->reco.pfp[id].shw.len);
  double this_endy = slc->reco.pfp[id].shw.start.y + (slc->reco.pfp[id].shw.dir.y * slc->reco.pfp[id].shw.len);
  double this_endz = slc->reco.pfp[id].shw.start.z + (slc->reco.pfp[id].shw.dir.z * slc->reco.pfp[id].shw.len);

  bool startx = ((-358.49 + 5. < slc->reco.pfp[id].shw.start.x) && (slc->reco.pfp[id].shw.start.x < -61.94 - 5.))
                || ((61.94 + 5. < slc->reco.pfp[id].shw.start.x) && (slc->reco.pfp[id].shw.start.x < 358.49 - 5. ));
  bool endx   = ((-358.49 + 5. < this_endx) && (this_endx < -61.94 - 5.))
                || ((61.94 + 5. < this_endx) && (this_endx < 358.49 - 5.));

  bool starty = (-181.86 + 5. < slc->reco.pfp[id].shw.start.y) && (slc->reco.pfp[id].shw.start.y < 134.96 - 5.);
  bool endy   = (-181.86 + 5. < this_endy) && (this_endy < 134.96 - 5.);

  bool startz = (-894.951 + 5. < slc->reco.pfp[id].shw.start.z) && (slc->reco.pfp[id].shw.start.z < 894.951 - 5.);
  bool endz   = (-894.951 + 5. < this_endz) && (this_endz < 894.951 - 5.);

  return (startx && endx && starty && endy && startz && endz);

}


////////////////////////
///////// Vars /////////
////////////////////////

const Var kTrkLen([](const caf::SRSliceProxy *slc) -> double {
  float len(-5.f);
  if (slc->reco.npfp > 0) {
    for(auto& pfp: slc -> reco.pfp) {
      const auto& trk = pfp.trk;
      if(!std::isnan(trk.len)) len = trk.len;
    }
  } 
  return len;
});

const MultiVar kTrkLenAll([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> len={};
  if (slc->reco.npfp > 0) {
    for(auto& pfp: slc -> reco.pfp) {
      const auto& trk = pfp.trk;
      if(!std::isnan(trk.len)) len.push_back(trk.len);
    }
  } 
  return len;
});

const Var kLeadingProtonId([](const caf::SRSliceProxy* slc) -> double {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  int muonIdx = kLeadingMuonId(slc);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    if ((int)i == muonIdx) continue;
    auto& p = slc->reco.pfp.at(i);
    if ( checkProton(p, slc) && p.trk.len > leading.second) {
      leading = std::make_pair(i, p.trk.len);
    }
  }
  // if limiting the proton energy
  /*if ( leading.first >= 0 ) {
    float leadingP = slc->reco.pfp.at(leading.first).trk.rangeP.p_proton;
    if ( leadingP < 0.4 || leadingP > 1.0) return -1; // leading proton momentum between [0.4, 1.0] GeV
  }*/
  return leading.first;
});

const Var kLeadingMuonId([](const caf::SRSliceProxy* slc) -> double {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    auto& p = slc->reco.pfp.at(i);
    if ( checkMuon(p, slc) && p.trk.len > leading.second) {
      leading = std::make_pair(i, p.trk.len);
    }
  }
  //cout << "returning leading muon " << leading.first << " out of " << slc->reco.pfp.size() << std::endl;
  return leading.first;
});

const Var kLargestShowerId([](const caf::SRSliceProxy* slc) -> double {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  int muonIdx = kLeadingMuonId(slc);
  int protonIdx = kLeadingProtonId(slc);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    if ((int)i == protonIdx || (int)i == muonIdx) continue;
    auto& p = slc->reco.pfp.at(i);
    if ( checkShower(p, slc) && p.shw.plane[2].energy > leading.second) {
      leading = std::make_pair(i, p.shw.plane[2].energy);
    }
  }
  return leading.first;
});

const Var kSecondLargestShowerId([](const caf::SRSliceProxy* slc) -> double {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  std::pair<int, float> secondary = std::make_pair(-1, -1);
  int muonIdx = kLeadingMuonId(slc);
  int protonIdx = kLeadingProtonId(slc);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    if ((int)i == protonIdx || (int)i == muonIdx) continue;
    auto& p = slc->reco.pfp.at(i);
    if ( checkShower(p, slc) && p.shw.plane[2].energy > leading.second) {
      secondary = std::make_pair(leading.first, leading.second);
      leading = std::make_pair(i, p.shw.plane[2].energy);
    }
  }
  return secondary.first;
});

const Var kLeadingChargedPionId([](const caf::SRSliceProxy* slc) -> double {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  int muonIdx = kLeadingMuonId(slc);
  int protonIdx = kLeadingProtonId(slc);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    if ((int)i == protonIdx || (int)i == muonIdx) continue;
    auto& p = slc->reco.pfp.at(i);
    if ( checkChargedPion(p, slc) && p.trk.len > leading.second) {
      leading = std::make_pair(i, p.trk.len);
    }
  }
  return leading.first;
});

const Var kLongestTrackId([](const caf::SRSliceProxy* slc) -> double {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    auto& p = slc->reco.pfp.at(i);
    if ( p.trk.len > leading.second) {
      leading = std::make_pair(i, p.trk.len);
    }
  }
  return leading.first;  
});

const Var kLongestTrack([](const caf::SRSliceProxy* slc) -> double {
  int id = kLongestTrackId(slc);
  return ( id >= 0 ? (double)slc->reco.pfp.at(id).trk.len : -5);
});

const Var kLeadingProtonMomentum([](const caf::SRSliceProxy *slc) -> double {
  //if (kLeadingProtonId(slc) >= 0) std::cout << "Proton momentum/id: " << slc->reco.pfp.at(kLeadingProtonId(slc)).trk.rangeP.p_proton << ", " << kLeadingProtonId(slc) << std::endl;
  return ( kLeadingProtonId(slc) >= 0 ? (double)slc->reco.pfp.at(kLeadingProtonId(slc)).trk.rangeP.p_proton : -5);
});

const MultiVar kProtonMomentum([](const caf::SRSliceProxy *slc) -> std::vector<double> { 
  std::vector<double> momenta;
  for (auto& p : slc->reco.pfp){
    if (checkProton(p, slc)) momenta.push_back(double(p.trk.rangeP.p_proton));
  }
  return momenta;
});

const Var kLeadingMuonMomentum([](const caf::SRSliceProxy *slc) -> double {
  return ( kLeadingMuonId(slc) >= 0 ? (double)slc->reco.pfp.at(kLeadingMuonId(slc)).trk.rangeP.p_muon : -5);
});

const MultiVar kMuonMomentum([](const caf::SRSliceProxy *slc) -> std::vector<double> { 
  std::vector<double> momenta;
  for (auto& p : slc->reco.pfp){
    if (checkMuon(p, slc)) momenta.push_back(double(p.trk.rangeP.p_muon));
  }
  return momenta;
});

const Var kTruth_NeutrinoE([](const caf::SRSliceProxy *slc) -> double {
  return ( kHasTruthMatch(slc) ? (float)slc->truth.E : -5.f );
});

const MultiVar kShowerE0([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 0) && double(s.shw.plane[0].energy > 0.02) ) e.push_back(double(s.shw.plane[0].energy));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar kShowerE1([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 1) && double(s.shw.plane[1].energy > 0.02) ) e.push_back(double(s.shw.plane[1].energy));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar kShowerE2([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 2) && double(s.shw.plane[2].energy > 0.02) ) e.push_back(double(s.shw.plane[2].energy));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar tShowerE_E0([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 0) && double(s.shw.plane[0].energy > 0.02) ) e.push_back(double(s.shw.truth.p.plane[0][0].visE));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar tShowerE_E1([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 1) && double(s.shw.plane[1].energy > 0.02) ) e.push_back(double(s.shw.truth.p.plane[0][1].visE));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar tShowerE_E2([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 2) && double(s.shw.plane[2].energy > 0.02) ) e.push_back(double(s.shw.truth.p.plane[0][2].visE));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar tShowerE_W0([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 0) && double(s.shw.plane[0].energy > 0.02) ) e.push_back(double(s.shw.truth.p.plane[1][0].visE));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar tShowerE_W1([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 1) && double(s.shw.plane[1].energy > 0.02) ) e.push_back(double(s.shw.truth.p.plane[1][1].visE));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const MultiVar tShowerE_W2([](const caf::SRSliceProxy *slc) -> std::vector<double> {
  std::vector<double> e = {};
  for (auto& s : slc->reco.pfp){
    if ( checkShower(s, slc, 2) && double(s.shw.plane[2].energy > 0.02) ) e.push_back(double(s.shw.truth.p.plane[1][2].visE));
  }
  if (e.size() == 0) e.push_back(-5.f); // need this or else Tree will break from non-equal entries
  return e;
});

const Var kLeadingShowerE([](const caf::SRSliceProxy *slc) -> double {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? (double)slc->reco.pfp.at(id).shw.plane[2].energy : -5.f );
});

const Var kSubLeadingShowerE([](const caf::SRSliceProxy *slc) -> double {
  int id = kSecondLargestShowerId(slc);
  return ( id >= 0 ? (double)slc->reco.pfp.at(id).shw.plane[2].energy : -5.f );
});

const Var kShowerAngle([](const caf::SRSliceProxy *slc) -> double {
  int primary = kLargestShowerId(slc);
  int second = kSecondLargestShowerId(slc);
  if (primary < 0 || second < 0) return -5;
  auto& trk1 = slc->reco.pfp.at(primary).shw;
  auto& trk2 =  slc->reco.pfp.at(second).shw;
  TVector3 dir1( trk1.dir.x, trk1.dir.y, trk1.dir.z );
  TVector3 dir2( trk2.dir.x, trk2.dir.y, trk2.dir.z );
	auto dir = dir1.Unit();
	double thGamma = dir.Angle(dir2.Unit());

  TVector3 vertex( slc->vertex.x, slc->vertex.y, slc->vertex.z );
  TVector3 start1( trk1.start.x, trk1.start.y, trk1.start.z );
  TVector3 start2( trk2.start.x, trk2.start.y, trk2.start.z );
  start1 = start1 - vertex;
  start2 = start2 - vertex;
	auto start = start1.Unit();
	double thGamma2 = start.Angle(start2.Unit());

  //std::cout << "checking angles: " << thGamma << ", " << thGamma2 << ", " << trk2.open_angle - trk1.open_angle << endl;

	auto cosTh = TMath::Cos(thGamma);  
  return (primary >= 0 && second >= 0 ? cosTh : -5);
});

const Var kShowerInvariantMass([](const caf::SRSliceProxy *slc) -> double {
  int primary = kLargestShowerId(slc);
  int second = kSecondLargestShowerId(slc);
  if (primary < 0 || second < 0) return -5;

  float e1 = slc->reco.pfp.at(primary).shw.plane[2].energy;
  float e2 = slc->reco.pfp.at(second).shw.plane[2].energy;
  //float e1 = slc->reco.pfp.at(primary).shw.truth.bestmatch.energy;
  //float e2 = slc->reco.pfp.at(second).shw.truth.bestmatch.energy;
  // optionally correct energy using fit function
  /*TF1* f = new TF1("f", "[0]*(TMath::Erf((x-[1])/([2]*sqrt(2))) - TMath::Erf((- [1])/([2]*sqrt(2))))/(1 - TMath::Erf((- [1])/([2]*sqrt(2))))", 0, 10);
  f->SetParameters(1.27596, 0.169703, 0.356733);
  e1 = f->Eval(e1);
  e2 = f->Eval(e2);*/

  auto& trk1 = slc->reco.pfp.at(primary).shw;
  auto& trk2 =  slc->reco.pfp.at(second).shw;

  TVector3 vertex( slc->vertex.x, slc->vertex.y, slc->vertex.z );
  TVector3 start1( trk1.start.x, trk1.start.y, trk1.start.z );
  TVector3 start2( trk2.start.x, trk2.start.y, trk2.start.z );
  // TVector3 start1( trk1.truth.p.start.x, trk1.truth.p.start.y, trk1.truth.p.start.z );
  // TVector3 start2( trk2.truth.p.start.x, trk2.truth.p.start.y, trk2.truth.p.start.z );
  start1 = start1 - vertex;
  start2 = start2 - vertex;
	auto start = start1.Unit();
	double thGamma = start.Angle(start2.Unit());
  auto cosTh = TMath::Cos(thGamma);

  float mass = std::sqrt(2*e1*e2*(1-cosTh)) * 1000.;
  
  //cout << TString(Form("Invariant mass: %f, energies %f, %f, angles %f, %f", mass, e1, e2, theta1, theta2)) << endl;
  //cout << TString(Form("\tInvariant mass: %f, energies %f, %f, angle %f", mass, e1, e2, thGamma)) << endl;

  return (primary >= 0 && second >= 0 ? mass : -5);
});

const Var kLeadingShowerDedx([](const caf::SRSliceProxy *slc) -> double {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? (double)slc->reco.pfp.at(id).shw.plane[2].dEdx : -5.f );
});

const Var kLeadingShowerConversionGap([](const caf::SRSliceProxy *slc) -> double  {
  int id = kLargestShowerId(slc);
  return ( id >=0 ? (double)slc->reco.pfp.at(id).shw.conversion_gap : -5.f );
});

const Var kLeadingShowerDensity([](const caf::SRSliceProxy *slc) -> double {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? (double)slc->reco.pfp.at(id).shw.density : -5.f );
});

const Var kLeadingShowerOpenAngle([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? (double)slc->reco.pfp.at(id).shw.open_angle : -5.f );
});

const Var kLeadingShowerCosAngle([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);  
  if ( id < 0 ) return -5.f;
  //double cosAngle = TMath::Cos(slc->reco.pfp.at(id).shw.open_angle);
  const auto& shwStart = slc->reco.pfp.at(id).shw.start;
  if (std::isnan(shwStart.x) || std::isinf(shwStart.x)) return -5.f;
  if (std::isnan(shwStart.y) || std::isinf(shwStart.y)) return -5.f;
  if (std::isnan(shwStart.z) || std::isinf(shwStart.z)) return -5.f;
  TVector3 showerPos = TVector3(shwStart.x, shwStart.y, shwStart.z);
  TVector3 zAxis = TVector3(0, 0, 1);
	double angle = zAxis.Angle(showerPos.Unit());
  double theta = TMath::Cos(angle);
  return ( id >= 0 ? theta : -5.f );
});

const Var kLeadingShowerEleP([](const caf::SRSliceProxy* slc) -> double { 
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -5.f;
  double e = (double)slc->reco.pfp.at(id).shw.plane[2].energy;
  double p = std::sqrt(std::pow(e, 2) - std::pow(e_mass, 2));
  return ( !std::isnan(p) ? p : -1.f );
});

const Var kNuScore([](const caf::SRSliceProxy *slc) -> double {
  return ( slc->nu_score );
});

// Emulated trigger time
const SpillVar kNuMISpillTriggerTime ( [](const caf::SRSpillProxy *sr) -> double {
  double triggerTime = 0.;

  double foundTriggerTime = sr->hdr.triggerinfo.trigger_within_gate;
  if ( !std::isnan(foundTriggerTime) && !std::isinf(foundTriggerTime) && foundTriggerTime < -15. ) triggerTime = -15.;
  else if ( std::isnan(foundTriggerTime) || std::isinf(foundTriggerTime) ) triggerTime = -16.;
  else if ( foundTriggerTime > 30. ) triggerTime = 30.; //1.6 for BNB
  else triggerTime = foundTriggerTime;

  return triggerTime;
});

const SpillVar kNTrueNue([](const caf::SRSpillProxy* spill) -> int {
    int goodNue = 0;
    for (const auto& nu : spill->mc.nu){
      bool isFV = false;
      if ( !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) ) {
        isFV = (( ( nu.position.x < -61.94 - 25 && nu.position.x > -358.49 + 25 ) ||
                 ( nu.position.x >  61.94 + 25 && nu.position.x <  358.49 - 25 )) &&
               ( ( nu.position.y > -181.86 + 25 && nu.position.y < 134.96 - 25 ) &&
                 ( nu.position.z > -894.95 + 30 && nu.position.z < 894.95 - 50 ) ));
      }
      if (std::abs(nu.initpdg) == 12 && isFV) goodNue++;
    }
    return goodNue;
});
  
  // Helper function
// int FindMaxShowerSliceIdx(const caf::SRSpillProxy* sr) {
//   int bestIdx = -1;
//   double maxEnergy = -1;
  
//   for (unsigned int i = 0; i < sr->slc.size(); i++) {
//     for (auto& pfp : sr->slc[i].reco.pfp) {
//       if (!checkShower(pfp, &sr->slc[i])) continue;
//       if (pfp.shw.plane[2].energy > maxEnergy) {
//         maxEnergy = pfp.shw.plane[2].energy;
//         bestIdx = i;
//       }
//     }
//   }
//   return bestIdx;
// }

// //SpillCut that caches the spill and computes max slice index
// const SpillCut kCacheMaxShowerSlice([](const caf::SRSpillProxy* sr) -> bool {
//   g_CurrentSpill = sr;
//   g_MaxShowerSliceIdx = FindMaxShowerSliceIdx(sr);
//   return (g_MaxShowerSliceIdx >= 0);  // Reject spills with no valid slice
// });

// //Slice Cut that checks if this is the max shower slice
// const Cut kIsMaxShowerEnergySlice([](const caf::SRSliceProxy* slc) -> bool {
//   if (!g_CurrentSpill || g_MaxShowerSliceIdx < 0) return false;
//   return (&g_CurrentSpill->slc[g_MaxShowerSliceIdx] == slc);
// });


const TruthVar kTGoodNue([](const caf::SRTrueInteractionProxy* nu) -> bool {
  bool isFV = false;
  if ( !std::isnan(nu->position.x) && !std::isnan(nu->position.y) && !std::isnan(nu->position.z) ) {
      isFV = (( ( nu->position.x < -61.94 - 25 && nu->position.x > -358.49 + 25 ) ||
                ( nu->position.x >  61.94 + 25 && nu->position.x <  358.49 - 25 )) &&
              ( ( nu->position.y > -181.86 + 25 && nu->position.y < 134.96 - 25 ) &&
                ( nu->position.z > -894.95 + 30 && nu->position.z < 894.95 - 50 ) ));
  }
  if (std::abs(nu->initpdg) == 12 && isFV) return true;
  return false;
});

const TruthVar kTGood1Mu1Pi0([](const caf::SRTrueInteractionProxy* nu) -> bool {
  bool isFV = false;
  if ( !std::isnan(nu->position.x) && !std::isnan(nu->position.y) && !std::isnan(nu->position.z) ) {
      isFV = (( ( nu->position.x < -61.94 - 25 && nu->position.x > -358.49 + 25 ) ||
                ( nu->position.x >  61.94 + 25 && nu->position.x <  358.49 - 25 )) &&
              ( ( nu->position.y > -181.86 + 25 && nu->position.y < 134.96 - 25 ) &&
                ( nu->position.z > -894.95 + 30 && nu->position.z < 894.95 - 50 ) ));
  }
  int muonCount = 0;
  int pionCount = 0;
  for(auto& p : nu->prim){
    if (std::abs(p.pdg) == 13) muonCount++;
    if (std::abs(p.pdg) == 111) pionCount++;
  }
  if (muonCount == 1 && pionCount == 1 && isFV) return true;
  return false;
});

const TruthVar kTnNue([](const caf::SRTrueInteractionProxy* nu) -> int {
  return std::isnan(nu->pdg) ? -5 : std::abs(nu->pdg) == 12;
});

const TruthVar kTnNumu([](const caf::SRTrueInteractionProxy* nu) -> int {
  return std::isnan(nu->pdg) ? -5 : std::abs(nu->pdg) == 14;
});

const Var kBaryDeltaZT([](const caf::SRSliceProxy* slc) -> double {
  return slc->barycenterFM.deltaZ_Trigger;
});

const Var kBaryRadiusT([](const caf::SRSliceProxy* slc) -> double {
  return slc->barycenterFM.radius_Trigger;
});

// vertex must be contined in detector
const Var kRFiducial([](const caf::SRSliceProxy* slc) -> bool {
  return containedVertex(slc->vertex, 25, 25, 25, 25, 50, 30);
});

// vertex must be contined in detector
const Var kTFiducial([](const caf::SRSliceProxy* slc) -> bool {
  return containedVertex(slc->truth.position, 25, 25, 25, 25, 50, 30);
});

const Var tGoodNue([](const caf::SRSliceProxy* slc) -> bool {
  if (slc->truth.index < 0) return false; // no truth match
  const auto& nu = slc->truth;
  if (std::isnan(nu.position.x) || std::isnan(nu.position.y) || std::isnan(nu.position.z)) return false; // check for nans

  bool isFV = false;
  if ( !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) ) {
      isFV = (( ( nu.position.x < -61.94 - 25 && nu.position.x > -358.49 + 25 ) ||
                ( nu.position.x >  61.94 + 25 && nu.position.x <  358.49 - 25 )) &&
              ( ( nu.position.y > -181.86 + 25 && nu.position.y < 134.96 - 25 ) &&
                ( nu.position.z > -894.95 + 30 && nu.position.z < 894.95 - 50 ) ));
  }
  if (std::abs(nu.initpdg) == 12 && isFV && !(nu.isnc)) return true;
  return false;
});

// shower contained 5cm from detector edges
const Var kShwContained([](const caf::SRSliceProxy* slc) {

  const int largestShwIdx(kLargestRecoShowerIdx(slc));
  if ( largestShwIdx==-1 ) return false;

  ///slc->rceo.pfp[largestShwIdx].shw
  double this_endx = slc->reco.pfp[largestShwIdx].shw.start.x + (slc->reco.pfp[largestShwIdx].shw.dir.x * slc->reco.pfp[largestShwIdx].shw.len);
  double this_endy = slc->reco.pfp[largestShwIdx].shw.start.y + (slc->reco.pfp[largestShwIdx].shw.dir.y * slc->reco.pfp[largestShwIdx].shw.len);
  double this_endz = slc->reco.pfp[largestShwIdx].shw.start.z + (slc->reco.pfp[largestShwIdx].shw.dir.z * slc->reco.pfp[largestShwIdx].shw.len);

  bool startx = ((-358.49 + 5. < slc->reco.pfp[largestShwIdx].shw.start.x) && (slc->reco.pfp[largestShwIdx].shw.start.x < -61.94 - 5.))
                || ((61.94 + 5. < slc->reco.pfp[largestShwIdx].shw.start.x) && (slc->reco.pfp[largestShwIdx].shw.start.x < 358.49 - 5. ));
  bool endx   = ((-358.49 + 5. < this_endx) && (this_endx < -61.94 - 5.))
                || ((61.94 + 5. < this_endx) && (this_endx < 358.49 - 5.));

  bool starty = (-181.86 + 5. < slc->reco.pfp[largestShwIdx].shw.start.y) && (slc->reco.pfp[largestShwIdx].shw.start.y < 134.96 - 5.);
  bool endy   = (-181.86 + 5. < this_endy) && (this_endy < 134.96 - 5.);

  bool startz = (-894.951 + 5. < slc->reco.pfp[largestShwIdx].shw.start.z) && (slc->reco.pfp[largestShwIdx].shw.start.z < 894.951 - 5.);
  bool endz   = (-894.951 + 5. < this_endz) && (this_endz < 894.951 - 5.);

  return (startx && endx && starty && endy && startz && endz);

});

const Var kMuonChi2([](const caf::SRSliceProxy* slc) -> double {
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? (double)slc->reco.pfp[id].trk.chi2pid[2].chi2_muon : -5);
});

const Var kProtonChi2([](const caf::SRSliceProxy* slc) -> double {
  int id = kLongestTrackId(slc);
  return (id >= 0 ? (double)slc->reco.pfp[id].trk.chi2pid[2].chi2_proton : -5);
});

const Var kNShowers([](const caf::SRSliceProxy* slc) -> int { 
  int nShowers = 0;
  for (auto& p : slc->reco.pfp){
    if ( checkShower(p, slc) ) nShowers++;
  }
  return nShowers;
});

const Var kMuonDeltaTrkVertex([](const caf::SRSliceProxy* slc) -> double {
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? trkDistance(slc->reco.pfp.at(id), slc) : -5 );
});

const Var kMuonCollectionHits([](const caf::SRSliceProxy* slc) -> int {
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? (int)slc->reco.pfp.at(id).trk.calo[2].nhit : -5 );
});

const Var kMuonTrkContained([](const caf::SRSliceProxy* slc) -> bool {
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? containedTrack(slc->reco.pfp.at(id), 5) : false);
});

const Var kInteractionMode([](const caf::SRSliceProxy* slc) -> int {
  return slc->truth.genie_mode;
});

const Var kNeutrinoType([](const caf::SRSliceProxy* slc) -> int {
  if ( slc->truth.index < 0 ) return 0; //cosmic or no truth match
  else if ( slc->truth.isnc==1 ) return 1; //NC
  else if ( abs(slc->truth.pdg) == 12 ) return 2; //nue
  else if ( abs(slc->truth.pdg) == 14 ) return 3; //numu
  else return -1;
});

const Var kShowerEResidual([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -5;
  double tEnergy = 0;
  if (  std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE ) && std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ) ) return -5;
  tEnergy = slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE > slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ? slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE : slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE;
  return ( slc->truth.index != -1 ? (slc->reco.pfp.at(id).shw.plane[2].energy - tEnergy) / slc->reco.pfp.at(id).shw.plane[2].energy : -5);
  
});

const Var kSubleadingShowerEResidual([](const caf::SRSliceProxy* slc) -> double {
  int id = kSecondLargestShowerId(slc);
  if ( id < 0 ) return -5;
  double tEnergy = 0;
  if (  std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE ) && std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ) ) return -5;
  tEnergy = slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE > slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ? slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE : slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE;
  return ( slc->truth.index != -1 ? (slc->reco.pfp.at(id).shw.plane[2].energy - tEnergy) / slc->reco.pfp.at(id).shw.plane[2].energy : -5);
  
});

const Var kTShowerE([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -5;
  double tEnergy = 0;
  if (  std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE ) && std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ) ) return -5;
  tEnergy = slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE > slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ? slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE : slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE;
  return ( slc->truth.index != -1 ? tEnergy : -5);
});

const Var kTSubleadingShowerE([](const caf::SRSliceProxy* slc) -> double {
  int id = kSecondLargestShowerId(slc);
  if ( id < 0 ) return -5;
  double tEnergy = 0;
  if (  std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE ) && std::isnan( slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ) ) return -5;
  tEnergy = slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE > slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE ? slc->reco.pfp.at(id).shw.truth.p.plane[0][2].visE : slc->reco.pfp.at(id).shw.truth.p.plane[1][2].visE;
  return ( slc->truth.index != -1 ? tEnergy : -5);
});

const Var kTShowerLength([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -5;
  double tLength = slc->reco.pfp.at(id).shw.truth.p.length;
  return ( slc->truth.index != -1 ? tLength : -5);
});

const Var kRShowerLength([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -5;
  return slc->reco.pfp.at(id).shw.len;
});

const Var kShowerLenResidual([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -5;
  double tLength = slc->reco.pfp.at(id).shw.truth.p.length;
  return ( slc->truth.index != -1 ? (slc->reco.pfp.at(id).shw.len - tLength) / slc->reco.pfp.at(id).shw.len : -5);
});

const Var kRShowerEndX([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -9999;
  double pt = slc->reco.pfp.at(id).shw.end.x;
  return ( std::isnan(pt) ? -9999 : pt );
});

const Var kRShowerEndY([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -9999;
  double pt = slc->reco.pfp.at(id).shw.end.y;
  return ( std::isnan(pt) ? -9999 : pt );
});

const Var kRShowerEndZ([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -9999;
  double pt = slc->reco.pfp.at(id).shw.end.z;
  return ( std::isnan(pt) ? -9999 : pt );
});

const Var kRShowerStartX([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -9999;
  double pt = slc->reco.pfp.at(id).shw.start.x;
  return ( std::isnan(pt) ? -9999 : pt );
});

const Var kRShowerStartY([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -9999;
  double pt = slc->reco.pfp.at(id).shw.start.y;
  return ( std::isnan(pt) ? -9999 : pt );
});

const Var kRShowerStartZ([](const caf::SRSliceProxy* slc) -> double {
  int id = kLargestShowerId(slc);
  if ( id < 0 ) return -9999;
  double pt = slc->reco.pfp.at(id).shw.start.z;
  return ( std::isnan(pt) ? -9999 : pt );
});

const Var kTrackLenResidual([](const caf::SRSliceProxy* slc) -> double {
  int id = kLongestTrackId(slc);
  if ( id < 0 ) return -9999;
  double tLength = slc->reco.pfp.at(id).trk.truth.p.length;
  return ( slc->truth.index != -1 ? (slc->reco.pfp.at(id).trk.len - tLength) / slc->reco.pfp.at(id).trk.len : -9999);
});

const Var kTTrackLen([](const caf::SRSliceProxy* slc) -> double {
  int id = kLongestTrackId(slc);
  if ( id < 0 ) return -5;
  double tLength = slc->reco.pfp.at(id).trk.truth.p.length;  
  return ( slc->truth.index != -1 ? tLength : -5 );
});

const Var kVertexResidual([](const caf::SRSliceProxy* slc) -> double {
  double dx = slc->vertex.x - slc->truth.position.x;
  double dy = slc->vertex.y - slc->truth.position.y;
  double dz = slc->vertex.z - slc->truth.position.z;
  double dVtx = std::sqrt(dx*dx + dy*dy + dz*dz);
  return (slc->truth.index != -1 ? dVtx : -9999);
});

// Get cryostat of vertex
const Var kVertexCryo([](const caf::SRSliceProxy* slc) -> int {
  return ( !std::isnan(slc->vertex.x) ? (slc->vertex.x < 0 ? 0 : 1) : -5.f );
});

const Var kTruthPDGID([](const caf::SRSliceProxy* slc) -> int {
  return ( slc->truth.index != -1 ? (int)slc->truth.pdg : -9999 );
});

const Var kCVNNueScore([](const caf::SRSliceProxy* slc) -> double {
  return ( !std::isnan(slc->cvn.nuescore) ? (float)slc->cvn.nuescore : -9999.f );
});

const Var kCVNNumuScore([](const caf::SRSliceProxy* slc) -> double {
  return ( !std::isnan(slc->cvn.numuscore) ? (float)slc->cvn.numuscore : -9999.f );
});

const Var kCVNNcScore([](const caf::SRSliceProxy* slc) -> double {
  return ( !std::isnan(slc->cvn.ncscore) ? (float)slc->cvn.ncscore : -9999.f );
});

const Var kCVNCosmicScore([](const caf::SRSliceProxy* slc) -> double {
  return ( !std::isnan(slc->cvn.cosmicscore) ? (float)slc->cvn.cosmicscore : -9999.f );
});

const Var rClearCosmic([](const caf::SRSliceProxy* slc) -> bool {
  return slc->is_clear_cosmic;
});

const Var rLargestShowerCut([](const caf::SRSliceProxy* slc) -> bool {
  return ( kLargestShowerId(slc) >= 0 );
});

////////////////////////
///////// Cuts /////////
////////////////////////

const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
  return !slc->is_clear_cosmic;
});

const Cut kIsClearCosmic([](const caf::SRSliceProxy* slc) {
  return slc->is_clear_cosmic;
});

const Cut kNucrlongtrkdiryCut([](const caf::SRSliceProxy* slc) {
  return slc->nuid.crlongtrkdiry > -0.7;
});

const Cut kLeadingProtonCut([](const caf::SRSliceProxy* slc) {
  return ( kLeadingProtonId(slc) >= 0 );
});

const Cut kLeadingMuonCut([](const caf::SRSliceProxy* slc) {
  return ( kLeadingMuonId(slc) >= 0 );
});

const Cut kLargestShowerCut([](const caf::SRSliceProxy* slc) {
  return ( kLargestShowerId(slc) >= 0 );
});

const Cut kLongestTrackCut([](const caf::SRSliceProxy* slc) {
  return ( kLongestTrackId(slc) >= 0 );
});

const Cut kGoodMuon([](const caf::SRSliceProxy* slc) {
  bool goodMuon = false;
  std::pair<int, float> leading = std::make_pair(-1, -1);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    auto& p = slc->reco.pfp.at(i);
    if ( checkMuon(p, slc) && p.trk.len > leading.second) {
      leading = std::make_pair(i, p.trk.len);
    }
  }
  if ( leading.first >= 0 ) return true;
  return false;
});

const Cut kGoodProton([](const caf::SRSliceProxy* slc) {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    auto& p = slc->reco.pfp.at(i);
    if ( checkProton(p, slc) && p.trk.len > leading.second) {
      leading = std::make_pair(i, p.trk.len);
    }
  }
  if ( leading.first >= 0 ) return true;
  return false;
});

const Cut kCVNNueCut([](const caf::SRSliceProxy* slc) {
  if (std::isnan(slc->cvn.nuescore)) return false;
  return ( (float)slc->cvn.nuescore > 0.85 );
});

/*const Cut kGoodEle([](const caf::SRSliceProxy* slc) {
  std::pair<int, float> leading = std::make_pair(-1, -1.);
  bool showerThreshold = false;
  for (std::size_t i(0); i < slc->reco.pfp.size(); ++i) {
    auto& p = slc->reco.pfp.at(i);
    if ( p.shw.bestplane_energy > 0.35 ) showerThreshold = true; // must be one shower with > 350 MeV
    break;
  }
  int largestShowerID = kLargestShowerId(slc);
  if ( largestShowerID >= 0 ){
    auto& p = slc->reco.pfp.at(kLargestShowerId(slc));
    if ( checkEle(p, slc) && showerThreshold ) return true;
  }
  return false;
});*/

// Fiducial Selections

// truth vertex contained in detector
const TruthCut kTFiducialCut([](const caf::SRTrueInteractionProxy* nu){
  return containedVertex(nu->position, 25, 25, 25, 25, 50, 30);
});

const TruthCut kTGoodNueCut([](const caf::SRTrueInteractionProxy* nu){
  return kTGoodNue(nu);
});

const TruthCut kTGood1Mu1Pi0Cut([](const caf::SRTrueInteractionProxy* nu) {
  return kTGood1Mu1Pi0(nu);
});

// vertex must be contined in detector
const Cut kRFiducialCut([](const caf::SRSliceProxy* slc) {
  return containedVertex(slc->vertex, 25, 25, 25, 25, 50, 30);
});

// longest track must be < 110 cm long
const Cut kLongestTrackLCut([](const caf::SRSliceProxy* slc) {
  int id = kLongestTrackId(slc);
  return (id >=0 ? slc->reco.pfp.at(id).trk.len < 190 : false);
});

// shower contained 5cm from detector edges
const Cut kShwContainedCut([](const caf::SRSliceProxy* slc) {

  const int id(kLargestShowerId(slc));
  return (id >= 0 ? ShowerContained(slc, id) : false);

});

const Cut kSecondaryShwContainedCut([](const caf::SRSliceProxy* slc) {
  int id = kSecondLargestShowerId(slc);
  return (id >= 0 ? ShowerContained(slc, id) : false);
});

  // Cut on having valid trigger time in approx. beam window
  const SpillCut kNuMIValidTrigger ( [](const caf::SRSpillProxy *sr) {
    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    return spillTriggerTime > -0.1 && spillTriggerTime < 9.7;
  });

  // Cut on having valid trigger time in approx. beam window
  const SpillCut kBNBValidTrigger ( [](const caf::SRSpillProxy *sr) {
    double spillTriggerTime = kNuMISpillTriggerTime(sr);
    return spillTriggerTime > -0.1 && spillTriggerTime < 1.6;
  });

// Truth Selections

  const SpillCut kTruth_VertexInFV([](const caf::SRSpillProxy* spill) {
    bool isFV = false;
    for (const auto& nu : spill->mc.nu){
      if ( !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) ) {
        isFV = (( ( nu.position.x < -61.94 - 25 && nu.position.x > -358.49 + 25 ) ||
                  ( nu.position.x >  61.94 + 25 && nu.position.x <  358.49 - 25 )) &&
                ( ( nu.position.y > -181.86 + 25 && nu.position.y < 134.96 - 25 ) &&
                  ( nu.position.z > -894.95 + 30 && nu.position.z < 894.95 - 50 ) ));
      }
    }
    return isFV;
  });

  const SpillCut kTruth_goodNue([](const caf::SRSpillProxy* spill) {
    bool isFV = false;
    bool goodNue = false;
    for (const auto& nu : spill->mc.nu){
      if ( !std::isnan(nu.position.x) && !std::isnan(nu.position.y) && !std::isnan(nu.position.z) ) {
        isFV = (( ( nu.position.x < -61.94 - 25 && nu.position.x > -358.49 + 25 ) ||
                 ( nu.position.x >  61.94 + 25 && nu.position.x <  358.49 - 25 )) &&
               ( ( nu.position.y > -181.86 + 25 && nu.position.y < 134.96 - 25 ) &&
                 ( nu.position.z > -894.95 + 30 && nu.position.z < 894.95 - 50 ) ));
      }
      if (std::abs(nu.initpdg) == 12 && isFV) goodNue = true;
    }
    return goodNue;
  });

// Cosmic Selections

const SpillCut kCRTPMT_correction([](const caf::SRSpillProxy* spill) {
	int nFlashesInBeam = 0;
	int nFlashesEntering = 0;

	for (const auto& match : spill->crtpmt_matches) {
		if (match.flashGateTime < -0.5 || match.flashGateTime > 9.6) continue;
			++nFlashesInBeam;
			bool hasEntering = false;
			for (const auto& crthit : match.matchedCRTHits) {
				if (crthit.PMTTimeDiff <= 0 && crthit.PMTTimeDiff > -1) {
					hasEntering = true;
					break;
				}
			}
		  if (hasEntering) ++nFlashesEntering;
	}
	if (nFlashesInBeam > 0 && nFlashesInBeam == nFlashesEntering) return false;
	return true;
});

const Cut kBaryDeltaZTCut([](const caf::SRSliceProxy* slc) {
    return slc->barycenterFM.deltaZ_Trigger >= 0. && slc->barycenterFM.deltaZ_Trigger < 100.;
  });

const Cut kBaryRadiusTCut([](const caf::SRSliceProxy* slc) {
  return slc->barycenterFM.radius_Trigger >= 0. && slc->barycenterFM.radius_Trigger < 100.;
});

// Electron Shower Selections

//largest shower energy > 200 MeV
const Cut kShowerEnergy([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? slc->reco.pfp.at(id).shw.plane[2].energy > 0.5 : false );
});

//largest shower dEdx < 3.625 MeV/cm
const Cut kShowerDedx([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? slc->reco.pfp.at(id).shw.plane[2].dEdx < 3.625 : false);
});

//largest shower conversion gap < 3.25 cm
const Cut kShowerConversionGap([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? slc->reco.pfp.at(id).shw.conversion_gap < 3.25 : false);
});

//largest shower density > 4.5 MeV/cm
const Cut kShowerDensity([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  return ( id >= 0 ? slc->reco.pfp.at(id).shw.density > 4.5 : false);
});

//cut on slice matched to a true neutrino interaction
const Cut kSliceNeutrinoMatched([](const caf::SRSliceProxy* slc) {
  return (slc->truth.index >= 0);
});

// 1mu1pi0X selections

bool t1Mu1Pi0(const caf::SRSliceProxy* slc, bool signal){
  // numu CC, 1 muon + 1 pi0 in primary interaction, fiducial volume
  if ( slc->truth.index < 0) return false;
  int goodMuons = 0;
  int goodPions = 0;
  for (const auto& p: slc->truth.prim){
    if ( std::abs(p.pdg) == 13 ) goodMuons++;
    if ( std::abs(p.pdg) == 111 ) goodPions++;
  }
  //cout << TString(Form("Good muons: %i, and good pi0: %i", goodMuons, goodPions)) << endl;
  //if (goodPions > 1) cout << "There are " << goodPions << " muons " << endl;
  if ( goodMuons != 1 ) return false;
  if ( goodPions < 1) return false;
  if ( signal && goodPions == 1 ) return true;
  if ( !signal && goodPions > 1 ) return true;
  return false;
}

const Cut kT1Mu1Pi0([](const caf::SRSliceProxy* slc) {
  //if (t1Mu1Pi0(slc, true)) cout << "good 1Mu1Pi0" << endl;
  return t1Mu1Pi0(slc, true);
});

const Cut kT1MuNPi0([](const caf::SRSliceProxy* slc) {
  //if (t1Mu1Pi0(slc, false)) cout << "good NMu1Pi0" << endl;
  return t1Mu1Pi0(slc, false);
});

// check muon candidate is primary with trk score > 0.45
const Cut kMuonTrk([](const caf::SRSliceProxy* slc) {
  int id = kLeadingMuonId(slc);
  if (id < 0) return false;
  return (slc->reco.pfp.at(id).trackScore > 0.45 && slc->reco.pfp.at(id).parent_is_primary);
});

const Cut kMuonDeltaTrkVertexCut([](const caf::SRSliceProxy* slc) {
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? trkDistance(slc->reco.pfp.at(id), slc) < 10 : false );
});

const Cut kMuonCollectionHitsCut([](const caf::SRSliceProxy* slc) {
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? slc->reco.pfp.at(id).trk.calo[2].nhit > 5 : false );
});

const Cut kMuonTrkContainedCut([](const caf::SRSliceProxy* slc) {
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? containedTrack(slc->reco.pfp.at(id), 5) : false);
});

const Cut kMuonChi2Cut([](const caf::SRSliceProxy* slc) {
  int id = kLeadingMuonId(slc);
  bool passChi2 = slc->reco.pfp.at(id).trk.chi2pid[2].chi2_muon < 30 && slc->reco.pfp.at(id).trk.chi2pid[2].chi2_proton > 60;
  bool contained = containedTrack(slc->reco.pfp.at(id), 5);
  if (id < 0) return false;
  return (contained ? passChi2 : true);
});

const Cut kMuonTrkLenCut([](const caf::SRSliceProxy* slc) { 
  int id = kLeadingMuonId(slc);
  return (id >= 0 ? slc->reco.pfp.at(id).trk.len > 50 : false);
});

const Cut goodPhotonCandidates([](const caf::SRSliceProxy* slc) {
  int primary = kLargestShowerId(slc);
  int secondary = kSecondLargestShowerId(slc);
  if (primary >= 0 && secondary >= 0) return true;
  return false;
});

const Cut primaryPhotonEnergyCut([](const caf::SRSliceProxy* slc) { 
  int id = kLargestShowerId(slc);
  return (id >= 0 ? slc->reco.pfp.at(id).shw.plane[2].energy > 0.03 : false);
});

const Cut secondaryPhotonEnergyCut([](const caf::SRSliceProxy* slc) { 
  int id = kSecondLargestShowerId(slc);
  return (id >= 0 ? slc->reco.pfp.at(id).shw.plane[2].energy > 0.03 : false);
});

const Cut nShowersCut([](const caf::SRSliceProxy* slc) { 
  int nShowers = 0;
  for (auto& p : slc->reco.pfp){
    if ( checkShower(p, slc) ) nShowers++;
  }
  //return (nShowers >= 2 && nShowers <= 3);
  return nShowers == 2;
});

const Cut openAngle([](const caf::SRSliceProxy* slc) { 
  int primary = kLargestShowerId(slc);
  int secondary = kSecondLargestShowerId(slc);
  if (primary < 0 || secondary < 0) return false;

  const auto& trk1 = slc->reco.pfp.at(primary).shw;
  const auto& trk2 = slc->reco.pfp.at(secondary).shw;

  TVector3 vertex( slc->vertex.x, slc->vertex.y, slc->vertex.z );
  TVector3 start1( trk1.start.x, trk1.start.y, trk1.start.z );
  TVector3 start2( trk2.start.x, trk2.start.y, trk2.start.z );
  start1 = start1 - vertex;
  start2 = start2 - vertex;
	auto start = start1.Unit();
	double thGamma = start.Angle(start2.Unit());
  return (thGamma > 0.35); //radians, corresponds to ~20 degrees (from microboone)
});

const Cut chargedPionVeto([](const caf::SRSliceProxy* slc) { 
  int leadingChargedPionIdx = kLeadingChargedPionId(slc);
  return (leadingChargedPionIdx >= 0 ? false : true);
});

const Cut kQECut([](const caf::SRSliceProxy* slc) {
  return (slc->truth.genie_mode == caf::kQE ? true : false);
});

const Cut kResCut([](const caf::SRSliceProxy* slc) {
  return (slc->truth.genie_mode == caf::kRes ? true : false);
});

const Cut kDisCut([](const caf::SRSliceProxy* slc) {
  return (slc->truth.genie_mode == caf::kDIS ? true : false);
});

const Cut kCohCut([](const caf::SRSliceProxy* slc) {
  return (slc->truth.genie_mode == caf::kCoh ? true : false);
});

const Cut kMECCut([](const caf::SRSliceProxy* slc) {
  return (slc->truth.genie_mode == caf::kMEC ? true : false);
});

const Cut kEastCryoCut([](const caf::SRSliceProxy* slc) {
  return slc->vertex.x < 0;
});

const Cut kWestCryoCut([](const caf::SRSliceProxy* slc) {
  return slc->vertex.x > 0;
});

const Cut kIsTrueElectron([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  if (id < 0) return false;
  return ( std::abs(slc->reco.pfp.at(id).shw.truth.p.pdg) == 11 );
});

const Cut kIsTrueMuon([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  if (id < 0) return false;
  return ( std::abs(slc->reco.pfp.at(id).shw.truth.p.pdg) == 13 );
});

const Cut kIsTruePhoton([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  if (id < 0) return false;
  return ( std::abs(slc->reco.pfp.at(id).shw.truth.p.pdg) == 22 );
});

const Cut kIsTrueProton([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  if (id < 0) return false;
  return ( std::abs(slc->reco.pfp.at(id).shw.truth.p.pdg) == 2212 );
});

const Cut kIsTrueNeutralPion([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  if (id < 0) return false;
  return ( slc->reco.pfp.at(id).shw.truth.p.pdg == 111 );
});

const Cut kIsTrueChargedPion([](const caf::SRSliceProxy* slc) {
  int id = kLargestShowerId(slc);
  if (id < 0) return false;
  return ( std::abs(slc->reco.pfp.at(id).shw.truth.p.pdg) == 211 );
});

const Cut kTrueNue([](const caf::SRSliceProxy* slc){
  return slc->truth.index >= 0 && slc->truth.pdg == 12 && slc->truth.isnc!=1;
});

const Cut kTrueANue([](const caf::SRSliceProxy* slc){
  return slc->truth.index >= 0 && slc->truth.pdg == -12 && slc->truth.isnc!=1;
});

const Cut kOtherModeCut = !kQECut && !kResCut && !kDisCut && !kCohCut && !kMECCut;
const Cut kOtherTruthPDG = !kIsTrueElectron && !kIsTrueMuon && !kIsTruePhoton && !kIsTrueProton && !kIsTrueNeutralPion && !kIsTrueChargedPion;

const Cut kShowersContained = kShwContainedCut && kSecondaryShwContainedCut;

const Cut kNueCC        = kIsNue && !kIsNC;
const Cut kNumuCC       = kIsNumu && !kIsNC;
const Cut kNC           = kIsNC;
const Cut kOther        = kIsNC;
const Cut kThisCosmic   = !kHasNu;

const Cut kBarycenterCuts = kBaryDeltaZTCut && kBaryRadiusTCut;

///////////////////////////////
///////// Binnings ///////////
//////////////////////////////
const Binning kSliceBinning     = Binning::Simple(3,0.f,3.f);

  
////////////////////////////
///////// Axes /////////////
////////////////////////////
const HistAxis axslice ("count",          kSliceBinning,          kCounting);




}


