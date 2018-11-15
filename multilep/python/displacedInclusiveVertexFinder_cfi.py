import FWCore.ParameterSet.Config as cms

displacedInclusiveVertexFinder  = cms.EDProducer("InclusiveVertexFinder",
       beamSpot = cms.InputTag("offlineBeamSpot"),
       primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
       #tracks = cms.InputTag("displacedAssocToTracks","displacedAssocToTracks","ANA"),
       tracks = cms.InputTag("unpackedTracksAndVertices"),
       minHits = cms.uint32(6), #old 8 -> 0 AOD produciton has problems with nhits
       maximumLongitudinalImpactParameter = cms.double(99999), #old  .3 -> infty
       minPt = cms.double(0.4), #old .8 -> 1 
       maxNTracks = cms.uint32(100), #old 30 -> 100

       clusterizer = cms.PSet(
           seedMax3DIPSignificance = cms.double(9999.),
           seedMax3DIPValue = cms.double(9999.),
           seedMin3DIPSignificance = cms.double(-9999), 
           seedMin3DIPValue = cms.double(-9999),
           clusterMaxDistance = cms.double(0.4), #500um #old .05 -> 1
           clusterMaxSignificance = cms.double(4.5), #4.5 sigma  #old  4.5 ---> infty
           distanceRatio = cms.double(20), # was cluster scale = 1 / density factor =0.05 
           clusterMinAngleCosine = cms.double(0.5), # only forward decays   #old accept backward decays (unboosted topologies) .5 -> -9999
       ),

       vertexMinAngleCosine = cms.double(0.95), # scalar prod direction of tracks and flight dir  #old accept backward decays  .95 -> .6 me 0
       vertexMinDLen2DSig = cms.double(2.5), #2.5 sigma 
       vertexMinDLenSig = cms.double(0.5), #0.5 sigma
       fitterSigmacut =  cms.double(3),
       fitterTini = cms.double(256),
       fitterRatio = cms.double(0.25),
       useDirectVertexFitter = cms.bool(True),
       useVertexReco  = cms.bool(True),
       vertexReco = cms.PSet(
               finder = cms.string('avr'),
               primcut = cms.double(1.0),
               seccut = cms.double(3),
               smoothing = cms.bool(True)
       )

)


