import uproot
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load ROOT file and TTree
file = uproot.open("../response.root")
tree = file["events"]

# Extract data as standard lists
E0_list = tree["E0"].array(library="np").tolist()
edep_list = tree["edep"].array(library="np").tolist()

# Filter where edep > 0
E0_filtered = []
edep_filtered = []

for e0, edeps in zip(E0_list, edep_list):
    for edep_val in edeps:
        if edep_val > 0:
            E0_filtered.append(e0)
            edep_filtered.append(edep_val)

# Define histogram parameters
xbins, xmin, xmax = 1500, 0, 150
ybins, ymin, ymax = 1500, 0, 150

# Create 2D histogram
H, xedges, yedges = np.histogram2d(
    E0_filtered, edep_filtered,
    bins=[xbins, ybins],
    range=[[xmin, xmax], [ymin, ymax]]
)

# Define a custom colormap where zero is white
cmap = plt.cm.viridis
cmap_with_white = cmap.copy()
cmap_with_white.set_bad(color='white')

# Mask zeros
H_masked = np.ma.masked_where(H.T == 0, H.T)

# Plot
plt.figure(figsize=(8, 6))
plt.imshow(
    H_masked,
    extent=[xmin, xmax, ymin, ymax],
    origin="lower",
    aspect="auto",
    cmap=cmap_with_white,
    interpolation='nearest'
)

plt.colorbar(label="Counts")
plt.xlabel("Photon energy (keV)")
plt.ylabel("Energy deposition (keV)")
plt.title("2D Histogram: Energy deposition (bin center) vs Photon energy (white = zero)")
plt.tight_layout()
plt.savefig("response_matrix.png", dpi=300)
#plt.show()

# Save histogram as CSV
xcenters = 0.5 * (xedges[1:] + xedges[:-1])
ycenters = 0.5 * (yedges[1:] + yedges[:-1])
df = pd.DataFrame(H.T, index=np.round(ycenters, 2), columns=np.round(xcenters, 2))
df.index.name = "Energy_Bin_Center "
df.columns.name = "Photon_Energy_Bin_Center"
df.to_excel("response.xlsx")

print("Plot saved as 'response.png' and CSV written to 'response.csv'")

