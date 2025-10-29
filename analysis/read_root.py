import uproot
import numpy as np

# Open ROOT file
file = uproot.open("2mev.root")

# Get the tree
tree = file["inp"]

# Get number of entries
nentries = tree.num_entries
print(f"Number of entries: {nentries}")

# Read branches as arrays (vectorized, no loop needed)
pos = tree["pos"].array()           # shape: (nentries, 3)
E0 = tree["E0"].array()
eventID = tree["eventID"].array()
v = tree["direction"].array()            #   direction, (px,py,pz)
energy = tree["energy"].array()
pdg = tree["pdg"].array()
parent = tree["parent"].array()

pdg_mapping={
        22:'gamma',
        11:'e-',
        -11:'e+'

        }

# Example: iterate over events if needed
for i in range(nentries):
    # Access i-th event
    pos_i = pos[i]
    E0_i = E0[i]
    energy_i = energy[i]
    print('====================')
    print(f'Hit #{i}')
    print("hit location:", pos[i]) #x,y,z
    print("hit energy:", energy_i)
    print("hit direction:", v[i])
    print("parent:", parent[i]) 
    # if parent ==0, it is a primary photon, otherwise secondary
    print("particle:", pdg_mapping.get(pdg[i],pdg[i]) )
    """
    -----------------
    particle,  code
    ------------------
    gamma, 22
    e- , 11
    e+, -11

    """



    if i>10:
        break

