import uproot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# -------------------------------
# User-defined parameters
# -------------------------------
bin_size_mm = 0.05  # 0.05 mm
# -------------------------------

# Open ROOT file and read tree
file = uproot.open("../response.root")
tree = file["inp"]

# Read x and y coordinates of hits
x = tree["pos"].array()[:, 0]  # pos[0]
y = tree["pos"].array()[:, 1]  # pos[1]
x=x.to_numpy()
y=y.to_numpy()
# Determine number of bins based on data range and bin size
x_min, x_max = -60, 60
y_min, y_max = -60,60
xbins = int((x_max - x_min) / bin_size_mm) + 1
ybins = int((y_max - y_min) / bin_size_mm) + 1

# Create 2D histogram
plt.figure(figsize=(8, 6))
plt.hist2d(
    x, y, bins=[xbins, ybins], 
    range=[[x_min, x_max], [y_min, y_max]], 
    cmap="viridis", 
# norm=LogNorm()  # <-- logarithmic color scale
)
plt.colorbar(label="Counts ")
plt.xlabel("X [mm]")
plt.ylabel("Y [mm]")
plt.title("2D Histogram of Hit Positions")
plt.tight_layout()
plt.savefig("hit_x_y_2d.png")
plt.show()

