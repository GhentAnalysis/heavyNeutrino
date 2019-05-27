#include "heavyNeutrino/multilep/interface/LeptonAnalyzer.h"

#include <map>


double reducedPdgId( int pdgId ){
    static const std::map< unsigned, double > pdgIdMap = {
        { 0, 0.},
        { 1, 0.125},
        { 2, 0.25},
        { 11, 0.375},
        { 13, 0.5},
        { 22, 0.625},
        { 130, 0.75},
        { 211, 0.875}
    };
    auto entry = pdgIdMap.find( fabs( pdgId ) );
    if( entry != pdgIdMap.cend() ){
        return entry->second;
    } else {
        return 1;
    }
}


double catchNanOrInf( double value ){
    if( std::isnan( value ) ){
        return -1;
    } else if( std::isinf( value ) ){
        return -1;
    } else{
        return value;
    }
}


void LeptonAnalyzer::fillClosestJetConstituents( const reco::Candidate& lepton, const pat::Jet* jet ){

    unsigned num_daughters = ( jet == nullptr ) ? 0 : jet->numberOfDaughters();
    if( num_daughters < maxJetSize ){
        for( unsigned i = num_daughters; i < maxJetSize; ++i ){
            _closestJetConstituentPt[_nL][i] = 0.;
            _closestJetConstituentEta[_nL][i] = 0.;
            _closestJetConstituentPhi[_nL][i] = 0.;
            _closestJetConstituentMass[_nL][i] = 0.;
            _closestJetConstituentPdgId[_nL][i] = 0;
            _closestJetConstituentPdgIdReduced[_nL][i] = 0.;
            _closestJetConstituentCharge[_nL][i] = 0;
            _closestJetConstituentdxy[_nL][i] = 0;
            _closestJetConstituentdxyError[_nL][i] = 0;
            _closestJetConstituentdxySig[_nL][i] = 0.;
            _closestJetConstituentdz[_nL][i] = 0;
            _closestJetConstituentdzError[_nL][i] = 0;
            _closestJetConstituentdzSig[_nL][i] = 0.;
            _closestJetConstituentNumberOfHits[_nL][i] = 0;
            _closestJetConstituentNumberOfPixelHits[_nL][i] = 0;
            _closestJetConstituentNumberOfPixelLayersWithMeasurement[_nL][i] = 0;
            _closestJetConstituentNumberOfStripLayersWithMeasurement[_nL][i] = 0;
            _closestJetConstituentHasTrack[_nL][i] = false;
            _closestJetConstituentVx[_nL][i] = 0.;
            _closestJetConstituentVy[_nL][i] = 0.;
            _closestJetConstituentVz[_nL][i] = 0.;
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
        _closestJetConstituentPdgIdReduced[_nL][d] = reducedPdgId( daughter->pdgId() );
        _closestJetConstituentCharge[_nL][d] = daughter->charge();
        _closestJetConstituentVx[_nL][d] = daughter->vx();
        _closestJetConstituentVy[_nL][d] = daughter->vy();
        _closestJetConstituentVz[_nL][d] = daughter->vz();
        if( daughter->hasTrackDetails() ){
            _closestJetConstituentdxy[_nL][d] = daughter->dxy();
            _closestJetConstituentdxyError[_nL][d] = catchNanOrInf( daughter->dxyError() );
            _closestJetConstituentdxySig[_nL][d] = catchNanOrInf( fabs( daughter->dxy()/daughter->dxyError() ) );
            _closestJetConstituentdz[_nL][d] = daughter->dz();
            _closestJetConstituentdzError[_nL][d] = catchNanOrInf( daughter->dzError() );
            _closestJetConstituentdzSig[_nL][d] = catchNanOrInf( fabs( daughter->dz()/daughter->dzError() ) );
            _closestJetConstituentNumberOfHits[_nL][d] = daughter->numberOfHits();
            _closestJetConstituentNumberOfPixelHits[_nL][d] = daughter->numberOfPixelHits();
            _closestJetConstituentNumberOfPixelLayersWithMeasurement[_nL][d] = daughter->pixelLayersWithMeasurement();
            _closestJetConstituentNumberOfStripLayersWithMeasurement[_nL][d] = daughter->stripLayersWithMeasurement();
            _closestJetConstituentHasTrack[_nL][d] = true;
        } else {
            _closestJetConstituentdxy[_nL][d] = -1;
            _closestJetConstituentdxyError[_nL][d] = 0;
            _closestJetConstituentdxySig[_nL][d] = -1.;
            _closestJetConstituentdz[_nL][d] = -1;
            _closestJetConstituentdzError[_nL][d] = 0;
            _closestJetConstituentdzSig[_nL][d] = -1.;
            _closestJetConstituentNumberOfHits[_nL][d] = -1;
            _closestJetConstituentNumberOfPixelHits[_nL][d] = -1;
            _closestJetConstituentNumberOfPixelLayersWithMeasurement[_nL][d] = -1;
            _closestJetConstituentNumberOfStripLayersWithMeasurement[_nL][d] = -1;
            _closestJetConstituentHasTrack[_nL][d] = false; 
        }
    }
}
