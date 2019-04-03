import FWCore.ParameterSet.Config as cms

#from HNL.DisplacedAdaptiveVertexFinder.unpackedTracksAndVertices_cfi import *
from heavyNeutrino.multilep.unpackedTracksAndVertices_cfi import *
from heavyNeutrino.multilep.displacedInclusiveVertexFinder_cfi import *
from heavyNeutrino.multilep.displacedVertexMerger_cfi import *
from heavyNeutrino.multilep.displacedTrackVertexArbitrator_cfi import *
#from heavyNeutrino.multilep.displacedSVAssociationIVF_cfi import * this module matches muons with vertices in Jessica and mohamed's code. It is not present in my code

displacedInclusiveSecondaryVertices = displacedVertexMerger.clone()
displacedInclusiveSecondaryVertices.secondaryVertices = cms.InputTag("displacedTrackVertexArbitrator")
displacedInclusiveSecondaryVertices.maxFraction = 0.05 #0.05 #old .2 -> .05
displacedInclusiveSecondaryVertices.minSignificance = 10


displacedInclusiveVertexing = cms.Sequence(unpackedTracksAndVertices * displacedInclusiveVertexFinder  * displacedVertexMerger * displacedTrackVertexArbitrator * displacedInclusiveSecondaryVertices ) 

