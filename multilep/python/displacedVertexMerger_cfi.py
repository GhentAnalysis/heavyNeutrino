import FWCore.ParameterSet.Config as cms

displacedVertexMerger = cms.EDProducer("VertexMerger",
       secondaryVertices = cms.InputTag("displacedInclusiveVertexFinder"),
       maxFraction = cms.double(0.7), #old .7 -> .5
       minSignificance = cms.double(2)) 


