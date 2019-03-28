#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"

void LeptonAnalyzer::fillClosestJetConstituents( const reco::Candidate& lepton, const pat::Jet* jet ){

    unsigned num_daughters = ( jet == nullptr ) ? 0 : jet->numberOfDaughters();
    if( num_daughters < maxJetSize ){
        for( unsigned i = num_daughters; i < maxJetSize; ++i ){
            _closestJetConstituentPt[i][_nL] = 0.;
            _closestJetConstituentEta[i][_nL] = 0.;
            _closestJetConstituentPhi[i][_nL] = 0.;
            _closestJetConstituentMass[i][_nL] = 0.;
            _closestJetConstituentPdgId[i][_nL] = 0;
            _closestJetConstituentCharge[i][_nL] = 0;
            _closestJetConstituentdxySig[i][_nL] = 0.;
            _closestJetConstituentdzSig[i][_nL] = 0.;
            _closestJetConstituentsNumberOfHits[i][_nL] = 0;
            _closestJetConstituentsNumberOfPixelHits[i][_nL] = 0;
            _closestJetConstituentsHasTrack[i][_nL] = false;
        }
    } 
    if( jet == nullptr ){
        _nClosestJetConstituents[_nL] = 0;
        return;
    }

    _nClosestJetConstituents[_nL] = jet->numberOfDaughters();
    for(unsigned d = 0; d < std::min( (unsigned) jet->numberOfDaughters(), maxJetSize); ++d){
        const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet->daughter(d);
        _closestJetConstituentPt[d][_nL] = daughter->pt();
        _closestJetConstituentEta[d][_nL] = daughter->eta();
        _closestJetConstituentPhi[d][_nL] = daughter->phi();
        _closestJetConstituentMass[d][_nL] = daughter->mass();
        _closestJetConstituentPdgId[d][_nL] = daughter->pdgId();

        _closestJetConstituentCharge[d][_nL] = daughter->charge();
        if( daughter->hasTrackDetails() ){
            _closestJetConstituentdxySig[d][_nL] = fabs( daughter->dxy()/daughter->dxyError() );
            _closestJetConstituentdzSig[d][_nL] = fabs( daughter->dz()/daughter->dzError() );
            _closestJetConstituentsNumberOfHits[d][_nL] = daughter->numberOfHits();
            _closestJetConstituentsNumberOfPixelHits[d][_nL] = daughter->numberOfPixelHits();
            _closestJetConstituentsHasTrack[d][_nL] = true;
        } else {
            _closestJetConstituentdxySig[d][_nL] = -1.;
            _closestJetConstituentdzSig[d][_nL] = -1.;
            _closestJetConstituentsNumberOfHits[d][_nL] = -1;
            _closestJetConstituentsNumberOfPixelHits[d][_nL] = -1;
            _closestJetConstituentsHasTrack[d][_nL] = false;
        }
    }
}
