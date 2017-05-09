#ifndef Tools_H
#define Tools_H

#include <vector>
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

namespace Tools{
	double getEffArea(const unsigned, const double);
	double getRelIso03(const pat::Muon&, const double);
	//double getIsoAlt(const pat::Muon&, double);
	double getRelIso03(const pat::Electron&, const double);
	double getRelIso(const reco::RecoCandidate&, const std::vector<pat::PackedCandidate>& , const double, const double);
	double getMiniIso(const reco::RecoCandidate&, const std::vector<pat::PackedCandidate>&, const double, const double);
}
#endif
