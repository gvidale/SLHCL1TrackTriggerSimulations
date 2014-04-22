import FWCore.ParameterSet.Config as cms

# Taken from https://github.com/gpetruc/cmssw/blob/micro-from700/PhysicsTools/PatAlgos/python/slimming/prunedGenParticles_cfi.py

prunedGenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
        "drop *", # this is the default
        "keep status == 3",  #keep event summary status3 (for pythia)
        "keep status == 21", #keep particles from the hard matrix element (for pythia8)
        "keep status == 22", #keep particles from the hard matrix element (for pythia8)
        "keep status == 23", #keep particles from the hard matrix element (for pythia8)
        "++keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15", # keep leptons and their ancestors
        "++keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16", # keep neutrinos and their ancestors
        "++keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25  || abs(pdgId) == 37", # keep W,Z,H and their ancestors
        "++keep pdgId == 22 && status == 1 && pt > 10",                    # keep gamma above 10 GeV
        "drop status == 2",                                                # drop the shower part of the history
        "keep++ abs(pdgId) == 15",                                         # but keep taus with their daughters
        #"++keep 4 <= abs(pdgId) <= 6 ",                                    # keep also heavy quarks
        #"++keep (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)",   # and their hadrons
        #"drop   status == 2 && abs(pdgId) == 21",                          # but remove again gluons in the inheritance chain
        "keep (4 <= abs(pdgId) <= 5) && (status = 2 || status = 11 || status = 71 || status = 72)", # keep heavy flavour quarks for parton-based jet flavour
        "keep (1 <= abs(pdgId) <= 3 || pdgId = 21) && (status = 2 || status = 11 || status = 71 || status = 72) && pt > 5", # keep light-flavour quarks and gluons for parton-based jet flavour
        "keep (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)", # keep b and c hadrons for hadron-based jet flavour
        "keep abs(pdgId) = 10411 || abs(pdgId) = 10421 || abs(pdgId) = 10413 || abs(pdgId) = 10423 || abs(pdgId) = 20413 || abs(pdgId) = 20423 || abs(pdgId) = 10431 || abs(pdgId) = 10433 || abs(pdgId) = 20433", # additional c hadrons for jet fragmentation studies
        "keep abs(pdgId) = 10511 || abs(pdgId) = 10521 || abs(pdgId) = 10513 || abs(pdgId) = 10523 || abs(pdgId) = 20513 || abs(pdgId) = 20523 || abs(pdgId) = 10531 || abs(pdgId) = 10533 || abs(pdgId) = 20533 || abs(pdgId) = 10541 || abs(pdgId) = 10543 || abs(pdgId) = 20543", # additional b hadrons for jet fragmentation studies
    )
)
