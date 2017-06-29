#!/usr/bin/env python

import json
from math import pi, sqrt, sin, sinh, asin, asinh, atan2
from itertools import izip
from collections import Counter
from ROOT import TFile, TTree

# Olmo Cerri version:
#   https://github.com/ocerri/SLHCL1TrackTriggerSimulations/blob/dev_tt25/TriggerTowerDefinition/script/TriggerTowerDefinitionChart.C
#   https://github.com/ocerri/SLHCL1TrackTriggerSimulations/blob/dev_tt25/TriggerTowerDefinition/script/TriggerTowerClass.h

# ______________________________________________________________________________
# Global objects

# Input single-muon file covering the full tracker
#stubs_input_file = "/uscms_data/d2/jiafu/L1TrackTrigger/CRAB_amsim_SLHC25p3_naper/tt48/stubs_tt48_50M.0.root"
# stubs_input_file = "root://cmsxrootd.fnal.gov//store/user/l1upgrades/SLHC/GEN/620_SLHC25p3_results/jftest1/ocerri_20161003/stubs_tt48_50M.0.root"
stubs_input_file = "root://cmsxrootd.fnal.gov//store/user/l1upgrades/SLHC/GEN/920/SingleMu/WholeTracker/ParticleGuns/Stubs/170613_044226/PGun_pt2_Tilted_TkOnly_Stubs.root"

# Configurations
nentries = 50e6
#nentries = 10e3
rStarPhi = 90.
rStarEta = 60.

# Map of TriggerTowers
trigger_towers = {}

# Map of ModuleConnections
module_connections = {}

# ______________________________________________________________________________
# Tiny functions

# Get layer number from moduleId
def decodeLayer(moduleId):
  return int(moduleId / 10000)

def deltaPhi(phi1, phi2):
  result = phi1 - phi2
  while result >   pi:  result -= 2*pi
  while result <= -pi:  result += 2*pi
  return result

def get_min_phi(phi1, phi2):
  if deltaPhi(phi1, phi2) < 0.:
    return phi1
  else:
    return phi2

def get_max_phi(phi1, phi2):
  if deltaPhi(phi1, phi2) < 0.:
    return phi2
  else:
    return phi1

# ______________________________________________________________________________
# Classes

class TriggerTowerModule:
  def __init__(self, moduleId):
    self.moduleId = moduleId
    self.connection = {}
    self.min_r   = 0.
    self.max_r   = 0.
    self.min_phi = 0.
    self.max_phi = 0.
    self.min_z   = 0.
    self.max_z   = 0.
    self.dirty   = True
  def add(self, tower, stub_r, stub_phi, stub_z):
    self.connection.setdefault(tower, 0)
    self.connection[tower] += 1
    if self.dirty:
      self.min_r   = stub_r
      self.max_r   = stub_r
      self.min_phi = stub_phi
      self.max_phi = stub_phi
      self.min_z   = stub_z
      self.max_z   = stub_z
      self.dirty   = False
    else:
      self.min_r   = min(self.min_r, stub_r)
      self.max_r   = max(self.max_r, stub_r)
      self.min_phi = get_min_phi(self.min_phi, stub_phi)
      self.max_phi = get_max_phi(self.max_phi, stub_phi)
      self.min_z   = min(self.min_z, stub_z)
      self.max_z   = max(self.max_z, stub_z)

class TriggerTowerLayer:
  def __init__(self, lay):
    self.lay     = lay
    self.min_r   = 0.
    self.max_r   = 0.
    self.min_phi = 0.
    self.max_phi = 0.
    self.min_z   = 0.
    self.max_z   = 0.
    self.dirty   = True
  def add(self, min_r, max_r, min_phi, max_phi, min_z, max_z):
    if self.dirty:
      self.min_r   = min_r
      self.max_r   = max_r
      self.min_phi = min_phi
      self.max_phi = max_phi
      self.min_z   = min_z
      self.max_z   = max_z
      self.dirty   = False
    else:
      self.min_r   = min(self.min_r, min_r)
      self.max_r   = max(self.max_r, max_r)
      self.min_phi = get_min_phi(self.min_phi, min_phi)
      self.max_phi = get_max_phi(self.max_phi, max_phi)
      self.min_z   = min(self.min_z, min_z)
      self.max_z   = max(self.max_z, max_z)

class TriggerTower:
  def __init__(self, tower):
    self.tower = tower
    self.modules = {}
    self.layers = {}
  def add_module(self, moduleId, stub_r, stub_phi, stub_z):
    self.modules.setdefault(moduleId, TriggerTowerModule(moduleId))
    self.modules[moduleId].add(self.tower, stub_r, stub_phi, stub_z)
  def find_boundaries(self):
    for moduleId, trigger_tower_module in self.modules.iteritems():
      lay = decodeLayer(moduleId)
      m = trigger_tower_module
      self.layers.setdefault(lay, TriggerTowerLayer(lay))
      self.layers[lay].add(m.min_r, m.max_r, m.min_phi, m.max_phi, m.min_z, m.max_z)

class ModuleConnection:
  def __init__(self, moduleId):
    self.moduleId = moduleId
    self.connection = {}
  def add_connection(self, connection):
    for tower, cnt in connection.iteritems():
      self.connection.setdefault(tower, 0)
      self.connection[tower] += cnt

# ______________________________________________________________________________
# Long functions

def get_phiStar_from_phi(phi, invPt, rStar=90.):
    mPtFactor = 0.3*3.8*1e-2/2.0
    if abs(mPtFactor * rStar * invPt) > 1.:
      raise Exception("track invPt=%f (pt=%f) cannot reach rStar=%f" % (invPt, abs(1.0/invPt), rStar))
    dphi = - asin(mPtFactor * rStar * invPt)
    phiStar = phi + dphi
    while phiStar < -pi:
      phiStar += pi*2.
    while phiStar >= pi:
      phiStar -= pi*2.
    return phiStar

def get_etaStar_from_eta(eta, z0, invPt, rStar=60.):
    mPtFactor = 0.3*3.8*1e-2/2.0
    if abs(rStar) < 1e-10: return eta
    if abs(invPt) < 1e-10 and invPt <  0.: invPt = -1e-10
    if abs(invPt) < 1e-10 and invPt >= 0.: invPt = +1e-10
    cot = sinh(eta)
    cotStar = (cot * (asin(mPtFactor * rStar * invPt)/(mPtFactor * invPt)) + z0) / rStar
    return asinh(cotStar)

def TrackParametersToTT(phi, invPt, eta, z0, apply_pt_cut=True):
    max_eta = 2.2
    max_z0 = 15.      # [cm]
    max_invPt = 1./3  # [1/GeV]
    etaStar = get_etaStar_from_eta(eta, z0, invPt)
    phiStar = get_phiStar_from_phi(phi, invPt)
    if abs(etaStar) > max_eta or abs(z0) > max_z0 or (apply_pt_cut and abs(invPt) > max_invPt):
      return -1

    tt_eta = int((etaStar + max_eta) / (max_eta*2./6))
    tt_phi = int((phiStar + pi) / (pi*2./8))
    tt_phi = (tt_phi + 4) % 8
    tt = tt_eta * 8 + tt_phi
    assert(0<=tt and tt<48)
    return tt

# Loop over input stubs.root, make trigger towers. Each trigger tower contains
# a list of moduleIds and the (r,phi,z) boundaries of the trigger tower
def make_trigger_towers():
  print "Making trigger towers ..."
  def load(tree, nentries=nentries):
    for ievt, evt in enumerate(tree):
      if (ievt == nentries):  break

      if (ievt % 1000 == 0):  print "Processing event: %i" % ievt

      pt     = evt.genParts_pt[0]
      eta    = evt.genParts_eta[0]
      phi    = evt.genParts_phi[0]
      vz     = evt.genParts_vz[0]
      charge = evt.genParts_charge[0]
      invPt  = float(charge) / pt

      tt = TrackParametersToTT(phi, invPt, eta, vz)
      if tt == -1:
        pass
      else:
        for (moduleId, stub_r, stub_phi, stub_z) in izip(evt.TTStubs_modId, evt.TTStubs_r, evt.TTStubs_phi, evt.TTStubs_z):
          trigger_towers.setdefault(tt, TriggerTower(tt))
          trigger_towers[tt].add_module(moduleId, stub_r, stub_phi, stub_z)
    return

  tfile = TFile.Open(stubs_input_file)
  tree = tfile.Get("ntupler/tree")
  load(tree)
  tfile.Close()
  return

# Make the connections of moduleIds to trigger towers. Each moduleId is
# connected to no more than 4 trigger towers
def make_module_connections():
  print "Making module connections ..."
  for tt, trigger_tower in trigger_towers.iteritems():
    for moduleId, trigger_tower_module in trigger_tower.modules.iteritems():
      module_connections.setdefault(moduleId, ModuleConnection(moduleId))
      module_connections[moduleId].add_connection(trigger_tower_module.connection)

  for moduleId, module_connection in module_connections.iteritems():
    most_common_4 = dict(Counter(module_connection.connection).most_common(4))
    module_connections[moduleId] = ModuleConnection(moduleId)
    module_connections[moduleId].add_connection(most_common_4)
  return

# Update the lists of moduleIds for all trigger towers based on the module
# connections
def update_trigger_towers():
  print "Updating trigger towers ..."
  for tt, trigger_tower in trigger_towers.iteritems():
    modules = {}
    for moduleId, trigger_tower_module in trigger_tower.modules.iteritems():
      connection = module_connections[moduleId].connection
      if tt in connection:
        modules[moduleId] = trigger_tower_module
    trigger_towers[tt].modules = modules
    trigger_towers[tt].find_boundaries()
  return

def write_trigger_sector_map():
  print "Writing trigger sector map ..."
  assert(len(trigger_towers) == 48)

  all_moduleIds = []
  writeme = []
  for tt in xrange(48):
    tt_moduleIds = [v.moduleId for k, v in trigger_towers[tt].modules.iteritems()]
    all_moduleIds += tt_moduleIds

    ieta = tt/8
    iphi = tt%8
    s = ",".join("{0}".format(n) for n in [ieta+1,iphi+1]+sorted(tt_moduleIds))
    #print s
    writeme.append(s)

  print "Num of trigger towers: %i" % len(trigger_towers)
  print "Num of unique modules: %i" % len(set(all_moduleIds))
  print "Num of all modules incl. duplicates: %i" % len(all_moduleIds)

  with open("trigger_sector_map_oc.csv", "w") as f:
    writeme = ["eta_idx, phi_idx, module_list"] + writeme
    writeme[-1] += "\n"  # add EOL
    f.write("\n".join(writeme))
  return

def write_trigger_sector_boundaries():
  print "Writing trigger sector boundaries ..."
  assert(len(trigger_towers) == 48)

  writeme = []
  for tt in xrange(48):
    sorted_iteritems = sorted(trigger_towers[tt].layers.iteritems(), key=lambda x: x[1].lay)
    for lay, trigger_tower_layer in sorted_iteritems:
      # Special attention to phi boundaries
      if trigger_tower_layer.max_phi < trigger_tower_layer.min_phi:
        trigger_tower_layer.max_phi += 2*pi

      tt_boundaries = [0] * 4
      if 5 <= lay <= 10:  # barrel
        tt_boundaries[0] = trigger_tower_layer.min_phi
        tt_boundaries[1] = trigger_tower_layer.max_phi
        tt_boundaries[2] = trigger_tower_layer.min_z
        tt_boundaries[3] = trigger_tower_layer.max_z
      elif 11 <= lay <= 15 or 18 <= lay <= 22:  # endcap
        tt_boundaries[0] = trigger_tower_layer.min_phi
        tt_boundaries[1] = trigger_tower_layer.max_phi
        tt_boundaries[2] = trigger_tower_layer.min_r
        tt_boundaries[3] = trigger_tower_layer.max_r

      alist = [tt,lay]+tt_boundaries
      s = "{0},{1},{2:.6f},{3:.6f},{4:.4f},{5:.4f}".format(*alist)
      #print s
      writeme.append(s)

  with open("trigger_sector_boundaries_oc.csv", "w") as f:
    writeme = ["tt/I,layer/I,phimin/D,phimax/D,zmin_cm/D,zmax_cm/D"] + writeme
    f.write("\n".join(writeme))
  return

def write_json_files():
  print "Writing json files ..."
  mymap = {}
  with open("trigger_sector_map_oc.csv", "r") as f:
    for line in f:
      if not line[0].isdigit():  continue
      values = line.split(",")

      # Convert to int
      values = [int(x) for x in values]

      key = (values[0]-1)*8 + (values[1]-1)
      values = sorted(values[2:])
      mymap[key] = values

  json.dump(mymap, open("trigger_sector_map_oc.json", "w"), sort_keys=True)

  with open("trigger_sector_boundaries_oc.csv", "r") as f:
    for line in f:
      if not line[0].isdigit():
        continue
      values = line.split(",")
      assert(len(values) == 6)

      # Convert to int or float
      values = [float(x) if "." in x else int(x) for x in values]

      key = values[0]*100 + values[1]
      values = values[2:]
      mymap[key] = values

  json.dump(mymap, open("trigger_sector_boundaries_oc.json", "w"), sort_keys=True)
  return

def generate_tt_definition():
  make_trigger_towers()
  make_module_connections()
  update_trigger_towers()
  write_trigger_sector_map()
  write_trigger_sector_boundaries()
  write_json_files()
  return

# ______________________________________________________________________________
# Main

if __name__ == '__main__':

  generate_tt_definition()
  print "DONE"
