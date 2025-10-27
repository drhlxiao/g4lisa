import uproot
import numpy as np

# Open ROOT file
file = uproot.open("response.root")

# Get the tree
tree = file["inp"]

# Get number of entries
nentries = tree.num_entries
print(f"Number of entries: {nentries}")

# Read branches as arrays (vectorized, no loop needed)
pos = tree["pos"].array()           # shape: (nentries, 3)
E0 = tree["E0"].array()
eventID = tree["eventID"].array()
itrack = tree["itrack"].array()
boundary = tree["boundary"].array()
pixelID = tree["pixelID"].array()
detectorID = tree["detectorID"].array()
v = tree["v"].array()               # shape: (nentries, 3)
theta = tree["theta"].array()
energy = tree["energy"].array()
pdg = tree["pdg"].array()
parent = tree["parent"].array()

# Example: iterate over events if needed
for i in range(nentries):
    # Access i-th event
    pos_i = pos[i]
    E0_i = E0[i]
    energy_i = energy[i]
    print("hit location:", pos[i])
    print("hit energy:", energy_i)
    # Process as needed
    # print(E0_i, pos_i, energy_i)

