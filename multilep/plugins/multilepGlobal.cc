#include "heavyNeutrino/multilep/plugins/multilepGlobal.h"
#include "heavyNeutrino/multilep/interface/Header.h"


multilepGlobal::multilepGlobal(const edm::ParameterSet& iConfig):
    vtxToken(                 consumes<std::vector<reco::Vertex>>(             iConfig.getParameter<edm::InputTag>("vertices"))),
    genEventInfoToken(        consumes<GenEventInfoProduct>(                   iConfig.getParameter<edm::InputTag>("genEventInfo"))),
    lheEventInfoToken(        consumes<LHEEventProduct>(                       iConfig.getParameter<edm::InputTag>("lheEventInfo"))),
    pileUpToken(              consumes<std::vector<PileupSummaryInfo>>(        iConfig.getParameter<edm::InputTag>("pileUpInfo"))),
    sampleIsData(                                                              iConfig.getUntrackedParameter<bool>("isData"))
{
    nVertices  = fs->make<TH1D>("nVertices", "Number of vertices", 120, 0, 120);

    if( sampleIsData ) return;

    nTrueInteractions   = fs->make<TH1D>("nTrueInteractions", "nTrueInteractions", 100, 0, 100);
    hCounter            = fs->make<TH1D>("hCounter",   "Events counter", 1, 0, 1);
    lheCounter          = fs->make<TH1D>("lheCounter", "Lhe weights",    110, 0, 110); //Counter to determine effect of pdf and scale uncertainties on the MC cross section
    psCounter           = fs->make<TH1D>("psCounter",  "Lhe weights",    14, 0, 14);
    tauCounter          = fs->make<TH1D>("tauCounter", "Number of taus", 3, 0, 3);
    HTCounter           = fs->make<TH1D>("HTCounter", "HT Incoming", 251, 0, 2510);
}


// ------------ method called once each job just before starting event loop  ------------
//void multilepGlobal::beginJob(edm::Service<TFileService>& fs){
//}


// ------------ method called for each event  ------------
void multilepGlobal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
    auto vertices       = getHandle(iEvent, vtxToken);
    unsigned _nVertex = vertices->size();
    if(_nVertex == 0) return;          // Don't consider 0 vertex events


    if( sampleIsData ){
        nVertices->Fill(_nVertex); 
    }else {

        edm::Handle<GenEventInfoProduct> genEventInfo          = getHandle(iEvent, genEventInfoToken);
        edm::Handle<LHEEventProduct> lheEventInfo              = getHandle(iEvent, lheEventInfoToken);
        edm::Handle<std::vector<PileupSummaryInfo>> pileUpInfo = getHandle(iEvent, pileUpToken);
        
        double _weight   = genEventInfo->weight();

        nVertices->Fill(_nVertex, _weight); 
        hCounter->Fill(0.5, _weight);
        nTrueInteractions->Fill(pileUpInfo->begin()->getTrueNumInteractions(), _weight);

        double _lheHTIncoming = 0.;
        unsigned _nTau = 0;

        if(!lheEventInfo.isValid()) return;


        // See http://home.thep.lu.se/~leif/LHEF/classLHEF_1_1HEPEUP.html for more detailes
        for(unsigned i = 0; i < (unsigned)lheEventInfo->hepeup().NUP; ++i){
            int _lheStatus          = lheEventInfo->hepeup().ISTUP[i];
            int _lhePdgId           = lheEventInfo->hepeup().IDUP[i];
            int _lheMother1         = lheEventInfo->hepeup().MOTHUP[i].first-1;
            int _lheMother2         = lheEventInfo->hepeup().MOTHUP[i].second-1;
            double px               = lheEventInfo->hepeup().PUP[i][0];
            double py               = lheEventInfo->hepeup().PUP[i][1];
            double pz               = lheEventInfo->hepeup().PUP[i][2];
            double _lheE            = lheEventInfo->hepeup().PUP[i][3];
            TLorentzVector vector   = TLorentzVector(px, py, pz, _lheE);
            double _lhePt           = vector.Et();

            bool hasIncomingMother = lheEventInfo->hepeup().ISTUP[_lheMother1]==-1 and lheEventInfo->hepeup().ISTUP[_lheMother2]==-1;  // Status -1 means mother is incoming
            bool quarkOrGluon      = (_lhePdgId==21 or (abs(_lhePdgId)>0 and abs(_lhePdgId) < 7));

            if(hasIncomingMother and _lheStatus==1 and quarkOrGluon) _lheHTIncoming += _lhePt; // To be used when an inclusive MC sample needs to be combined with HT-splitted samples
            if(abs(_lhePdgId)==15) ++_nTau;
        }

        tauCounter->Fill(_nTau, _weight);
        HTCounter->Fill(_lheHTIncoming, _weight);

        //Store LHE weights to compute pdf and scale uncertainties, as described on https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
        for(unsigned i = 0; i < std::min( (unsigned) 110, (unsigned) lheEventInfo->weights().size()); ++i){// 110 for MC@NLO, 254 for powheg, 446(!) for madgraph, 0 for some old samples,...
            double _lheWeight = lheEventInfo->weights()[i].wgt/lheEventInfo->originalXWGTUP();
            lheCounter->Fill(i + 0.5, _lheWeight*_weight);
        }

        //some tests for PS weight extraction
        const std::vector<double>& psWeights = genEventInfo->weights();
        unsigned _nPsWeights = std::min( (unsigned) 14, (unsigned) psWeights.size() );
        for(unsigned ps = 0; ps < _nPsWeights; ++ps){
            double _psWeight = psWeights[ps]/_weight;
            psCounter->Fill(ps + 0.5, _psWeight*_weight);
        }
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(multilepGlobal);
