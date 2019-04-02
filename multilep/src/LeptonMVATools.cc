#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"

void LeptonAnalyzer::fillClosestJetConstituents( const reco::Candidate& lepton, const pat::Jet* jet ){

    unsigned num_daughters = ( jet == nullptr ) ? 0 : jet->numberOfDaughters();
    if( num_daughters < maxJetSize ){
        for( unsigned i = num_daughters; i < maxJetSize; ++i ){
            _closestJetConstituentPt[_nL][i] = 0.;
            _closestJetConstituentEta[_nL][i] = 0.;
            _closestJetConstituentPhi[_nL][i] = 0.;
            _closestJetConstituentMass[_nL][i] = 0.;
            _closestJetConstituentPdgId[_nL][i] = 0;
            _closestJetConstituentCharge[_nL][i] = 0;
            _closestJetConstituentdxySig[_nL][i] = 0.;
            _closestJetConstituentdzSig[_nL][i] = 0.;
            _closestJetConstituentsNumberOfHits[_nL][i] = 0;
            _closestJetConstituentsNumberOfPixelHits[_nL][i] = 0;
            _closestJetConstituentsHasTrack[_nL][i] = false;
        }
    } 
    if( jet == nullptr ){
        _nClosestJetConstituents[_nL] = 0;
        return;
    }

    _nClosestJetConstituents[_nL] = jet->numberOfDaughters();
    for(unsigned d = 0; d < std::min( (unsigned) jet->numberOfDaughters(), maxJetSize); ++d){
        const pat::PackedCandidate* daughter = (const pat::PackedCandidate*) jet->daughter(d);
        _closestJetConstituentPt[_nL][d] = daughter->pt();
        _closestJetConstituentEta[_nL][d] = daughter->eta();
        _closestJetConstituentPhi[_nL][d] = daughter->phi();
        _closestJetConstituentMass[_nL][d] = daughter->mass();
        _closestJetConstituentPdgId[_nL][d] = daughter->pdgId();

        _closestJetConstituentCharge[_nL][d] = daughter->charge();
        if( daughter->hasTrackDetails() ){
            _closestJetConstituentdxySig[_nL][d] = fabs( daughter->dxy()/daughter->dxyError() );
            _closestJetConstituentdzSig[_nL][d] = fabs( daughter->dz()/daughter->dzError() ); 
            _closestJetConstituentsNumberOfHits[_nL][d] = daughter->numberOfHits();
            _closestJetConstituentsNumberOfPixelHits[_nL][d] = daughter->numberOfPixelHits();
            _closestJetConstituentsHasTrack[_nL][d] = true;
        } else {
            _closestJetConstituentdxySig[_nL][d] = -1.;
            _closestJetConstituentdzSig[_nL][d] = -1.;
            _closestJetConstituentsNumberOfHits[_nL][d] = -1;
            _closestJetConstituentsNumberOfPixelHits[_nL][d] = -1;
            _closestJetConstituentsHasTrack[_nL][d] = false; 
        }
    }
}
