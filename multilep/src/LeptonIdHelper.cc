#include "heavyNeutrino/multilep/interface/LeptonIdHelper.h"

LeptonIdHelper::LeptonIdHelper(const LeptonAnalyzer& lepAn, const unsigned index, const unsigned flav){
    flavor = flav;
    pt = lepAn._lPt[index];
    eta = lepAn._lEta[index];
    miniIsoCharged = lepAn._miniIsoCharged[index];
    miniIso = lepAn._miniIso[index];
    relIso = lepAn._relIso[index];
    ptRel = lepAn._ptRel[index];
    ptRatio = lepAn._ptRatio[index];
    closestJetCsvV2 = lepAn._closestJetCsvV2[index];
    sip3d = lepAn._3dIPSig[index];
    dxy = lepAn._dxy[index];
    dz = lepAn._dz[index];
    if(flavor == 0){
        eleMva = lepAn._lElectronMva[index];
        eleTrigEmu = lepAn._lElectronPassEmu[index];
    } else if(flavor == 1){
        mediumMuon = lepAn._lPOGMedium[index];              //The medium muon variable in leptonAnalyzer might have to be renamed
    } else if(flavor == 2){
        ;                                                   //Tau id to be added later
    }
}
