
import sys
import uproot
import numpy as np
import matplotlib.pyplot as plt


def process_file(fname):
    # Open the ROOT file and access the tree
    file = uproot.open(fname)
    tree = file["inp"]

    # Read branches into arrays
    pos = tree["pos"].array(library="np")
    E0 = tree["E0"].array(library="np")
    eventID = tree["eventID"].array(library="np")
    direction = tree["direction"].array(library="np")
    energy = tree["energy"].array(library="np")
    pdg = tree["pdg"].array(library="np")
    parent = tree["parent"].array(library="np")
    #theta= tree["theta"].array(library="np")

    dx = direction[:, 0]
    dy = direction[:, 1]
    dz = direction[:, 2]
    theta=np.rad2deg( np.atan2(np.sqrt(dx**2 + dy**2), dz))


    # Extract x, y, z coordinates from pos array
    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]

    print(f"Total number of entries: {len(energy)}")
    print(f"Primary photons (parent=0): {np.sum(parent == 0)}")
    print(f"Secondary particles (parent>0): {np.sum(parent > 0)}")

    # Create figure with subplots
    fig = plt.figure(figsize=(15, 10))

    # 1) 2D histogram showing hit location distribution (x-y plane)
    ax1 = fig.add_subplot(2, 2, 1)
    h1 = ax1.hist2d(x, y, bins=1000, cmap='viridis')
    ax1.set_xlabel('X Position (mm)')
    ax1.set_ylabel('Y Position (mm)')
    ax1.set_title('Hit Location Distribution (X-Y plane)')
    plt.colorbar(h1[3], ax=ax1, label='Counts')


    # 2) Energy spectrum of all hits
    ax3 = fig.add_subplot(2, 2, 2)
    ax3.hist(energy, bins=100, color='blue', alpha=0.7, edgecolor='black')
    ax3.set_xlabel('Energy (MeV)')
    ax3.set_ylabel('Counts')
    ax3.set_title('Energy Spectrum of All Hits')
    ax3.set_yscale('log')
    ax3.grid(True, alpha=0.3)


    # --- Mapping ---
    pdg_mapping = {
        22: 'gamma',
        11: 'e-',
        -11: 'e+'
    }

    particle_type_mapping = {
        22: 0,
        11: 1,
        -11: 2
    }

    # Convert pdg array to integer codes
    to_type = np.array([particle_type_mapping.get(p, 3) for p in pdg])

    ax4 = fig.add_subplot(2, 2, 3)
    pdg_counts = np.bincount(to_type, minlength=4)  # ensure enough bins for code=3 (unknown)
    ax4.bar(range(len(pdg_counts)), pdg_counts, color='green', alpha=0.7, edgecolor='black')
    ax4.set_xlabel('Particle Type')
    ax4.set_ylabel('Counts')
    ax4.set_yscale('log')
    # Set tick labels
    tick_labels = [pdg_mapping.get(k, 'other') for k in [22, 11, -11, 0]]  # 0 is placeholder for unknown
    ax4.set_xticks(range(len(pdg_counts)))
    ax4.set_xticklabels(tick_labels)

    ax4.grid(True, alpha=0.3)

    # 4) Energy spectrum of secondaries (parent > 0)
    theta_mask = (pdg==22)
    ax5 = fig.add_subplot(2, 2, 4)
    ax5.hist(theta[theta_mask], bins=100, color='red', alpha=0.7, edgecolor='black')
    ax5.set_xlabel('Theta (deg)')
    ax5.set_ylabel('Counts')
    ax5.set_yscale('log')
    ax5.set_title(f'Distribution of Incident Angle at the detector')
    ax5.grid(True, alpha=0.3)



    plt.tight_layout()
    plt.savefig('analysis_results.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__=='__main__':
    if len(sys.argv)!=2:
        print('usage: python analyze.py <filename>')
        sys.exit(1)
    process_file(sys.argv[1])
