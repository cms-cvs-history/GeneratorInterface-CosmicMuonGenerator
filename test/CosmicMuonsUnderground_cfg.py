import FWCore.ParameterSet.Config as cms

process = cms.Process("runCosMuoGen")
process.load("GeneratorInterface.CosmicMuonGenerator.CMSCGENsource_cfi")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    sourceSeed = cms.untracked.uint32(135799468)
)

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(500)
    input = cms.untracked.int32(20)
#    input = cms.untracked.int32(100000)
)
process.CMSCGEN_out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('cosmic.root')
)

process.outpath = cms.EndPath(process.CMSCGEN_out)

#process.CosMuoGenSource.MinP = 10.
process.CosMuoGenSource.MinP = 3.
#process.CosMuoGenSource.MinP = 300.

#process.CosMuoGenSource.MaxTheta = 80.
process.CosMuoGenSource.MaxTheta = 89.

# Plug z-position [mm] (default=-14000. = on Shaft) 
#process.CosMuoGenSource.PlugVz = 5000.;
#process.CosMuoGenSource.PlugVz = -33000.;

# z-position of centre of target cylinder [mm] (default=0.)
#process.CosMuoGenSource.ZCentrOfTarget = 0.
#process.CosMuoGenSource.ZCentrOfTarget = 15000. 
#process.CosMuoGenSource.TrackerOnly = True

#Read in Multi muon events or generate single muon events (MultiMuon=false = default)
process.CosMuoGenSource.MultiMuon = True
#process.CosMuoGenSource.MultiMuonNmin = 1
#process.CosMuoGenSource.MultiMuonNmin = 2
process.CosMuoGenSource.MultiMuonNmin = 3
#process.CosMuoGenSource.MultiMuonFileName = "MultiEventsIn.root"
#process.CosMuoGenSource.MultiMuonFileName = "test_25gev.root"
#process.CosMuoGenSource.MultiMuonFileName = "test_150gev.root"
#process.CosMuoGenSource.MultiMuonFileName = "test_800gev.root"
#process.CosMuoGenSource.MultiMuonFileName = "test_3000gev.root"
##process.CosMuoGenSource.MultiMuonFileName = "CORSIKA6900_800_3000GeV_1k.root"
####process.CosMuoGenSource.MultiMuonFileName = "CORSIKA6900_800_3000GeV_10k.root"
process.CosMuoGenSource.MultiMuonFileName = "CORSIKA6900_3_10TeV_100k.root"
#process.CosMuoGenSource.MultiMuonFileName = "CORSIKA6900_800_3000GeV_100k.root"
###process.CosMuoGenSource.MultiMuonFileFirstEvent = 1
 


#process.CosMuoGenSource.Verbosity = True
