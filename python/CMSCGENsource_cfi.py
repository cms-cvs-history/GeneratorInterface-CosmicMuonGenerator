import FWCore.ParameterSet.Config as cms

source = cms.Source("CosMuoGenSource",
    RadiusOfTarget = cms.double(8000.0),
    TIFOnly_constant = cms.bool(False),
    TIFOnly_linear = cms.bool(False),
    ElossScaleFactor = cms.double(1.0),
    MaxPhi = cms.double(360.0),
    MaxTheta = cms.double(84.26),
    # MaxTheta = cms.double(89.0),
    MaxT0 = cms.double(12.5),
    ZDistOfTarget = cms.double(15000.0),
    ZCentrOfTarget = cms.double(0.0),
    MinP = cms.double(10.0),
    MinT0 = cms.double(-12.5),
    MTCCHalf = cms.bool(False),
    TrackerOnly = cms.bool(False),
    MultiMuon = cms.bool(False),
    #MultiMuon = cms.bool(True),
    MultiMuonFileName = cms.string("test_3000gev.root"),
    MultiMuonFileFirstEvent = cms.int32(1),
    MultiMuonNmin = cms.int32(2),                
    MinTheta = cms.double(0.0),
    MinP_CMS = cms.double(-1.0), ##negative means MinP_CMS = MinP. Only change this if you know what you are doing!

    MaxP = cms.double(3000.0),
    MinPhi = cms.double(0.0),
    PlugVx = cms.double(0.0),
    PlugVz = cms.double(-14000.0),                
    Verbosity = cms.bool(False)
)


 
