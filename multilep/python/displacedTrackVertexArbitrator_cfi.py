import FWCore.ParameterSet.Config as cms

displacedTrackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
       beamSpot = cms.InputTag("offlineBeamSpot"),
       primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
       #tracks = cms.InputTag("displacedAssocToTracks","displacedAssocToTracks","ANA"),
       tracks = cms.InputTag("unpackedTracksAndVertices"),
       secondaryVertices = cms.InputTag("displacedVertexMerger"),
       dLenFraction = cms.double(0.333), #old .333 -> .2
       dRCut = cms.double(1), # old .4 -> 1   me 3
       distCut = cms.double(0.1), #old .04 -> .1 
       sigCut = cms.double(5), #old 5->10 
       fitterSigmacut =  cms.double(3),
       fitterTini = cms.double(256),
       fitterRatio = cms.double(0.25),
       trackMinLayers = cms.int32(4), #old 4-> 0
       trackMinPt = cms.double(.4), 
       trackMinPixels = cms.int32(0) #old 1 -> 0
)


