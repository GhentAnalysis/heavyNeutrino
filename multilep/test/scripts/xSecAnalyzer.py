#
# Cross section analyzer, use as
# cmsRun xSecAnalyzer.py inputFiles="file:///pnfs/iihe/cms/.../file"
#
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.parseArguments()

process = cms.Process('xSectionAnalyzer')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source(
    "PoolSource",
    fileNames  = cms.untracked.vstring(options.inputFiles),
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
)


process.dummy = cms.EDAnalyzer("GenXSecAnalyzer")

process.ana = cms.Path(process.dummy)
process.schedule = cms.Schedule(process.ana)
