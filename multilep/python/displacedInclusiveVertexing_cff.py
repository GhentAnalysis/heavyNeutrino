import FWCore.ParameterSet.Config as cms

#from HNL.DisplacedAdaptiveVertexFinder.unpackedTracksAndVertices_cfi import *
from PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi import *
from HNL.DisplacedAdaptiveVertexFinder.displacedInclusiveVertexFinder_cfi import *
from HNL.DisplacedAdaptiveVertexFinder.displacedVertexMerger_cfi import *
from HNL.DisplacedAdaptiveVertexFinder.displacedTrackVertexArbitrator_cfi import *
from HNL.DisplacedSVAssociator.displacedSVAssociationIVF_cfi import *

displacedInclusiveSecondaryVertices = displacedVertexMerger.clone()
displacedInclusiveSecondaryVertices.secondaryVertices = cms.InputTag("displacedTrackVertexArbitrator")
displacedInclusiveSecondaryVertices.maxFraction = 0.05 #0.05 #old .2 -> .05
displacedInclusiveSecondaryVertices.minSignificance = 10


displacedInclusiveVertexing = cms.Sequence(unpackedTracksAndVertices * displacedInclusiveVertexFinder  * displacedVertexMerger * displacedTrackVertexArbitrator * displacedInclusiveSecondaryVertices ) 

