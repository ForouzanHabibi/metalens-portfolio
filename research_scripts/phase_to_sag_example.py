import numpy as np
import matplotlib.pyplot as plt
import os
import importlib.util
import sys
from skimage.restoration import unwrap_phase  # ✅ Required for phase unwrapping
import scipy.io as sio


# Define the path to the lumapi.py file
lumapi_path = r'C:\Program Files\Lumerical\v241\api\python'
spec = importlib.util.spec_from_file_location("lumapi", os.path.join(lumapi_path, "lumapi.py"))
lumapi = importlib.util.module_from_spec(spec)
sys.modules["lumapi"] = lumapi
spec.loader.exec_module(lumapi)
import lumapi  # Lumerical API
import scipy.io as sio
plt.rcdefaults()  # Reset all matplotlib settings to default
plt.set_cmap('jet')
# Path to your Lumerical file
lumfilename = "metalens_from_p1000_d_f_T_phase_R_A_LCP.txt_H_on_RCP_R15.0um_f10.0um_Marco approach_2025-02-24_17-01-52.fsp" # RCP light on LCP LUT

# Define simulation parameters
wavelength = 1.55e-6  # meters
n_material = 2.2       # Refractive index of BIG material
R = 15e-6              # Lens radius in meters
height = 1300e-9       # Thickness of metalens

# Open Lumerical FDTD session
with lumapi.FDTD() as fdtd:
    fdtd.load(lumfilename)

    # Extract near-field data from DFT monitor at output plane
    Enf = fdtd.getresult("DFTmonitor0", "E")
    Hnf = fdtd.getresult("DFTmonitor0", "H")

    x = Enf['x'][:, 0]
    y = Enf['y'][:, 0]

    # Extract Ex, Ey, Ez components
    Enfx = Enf['E'][:, :, 0, :, 0]
    Enfy = Enf['E'][:, :, 0, :, 1]
    Enfz = Enf['E'][:, :, 0, :, 2]

    # Combine fields
    Enf_xy = Enfx[:, :, 0] + Enfy[:, :, 0] + Enfz[:, :, 0]

    # Calculate raw phase and unwrap
    raw_phase = np.angle(Enf_xy)
    unwrapped_phase = unwrap_phase(raw_phase)

    # Convert phase to surface sag for Zemax
    sag_map = (unwrapped_phase * wavelength) / (2 * np.pi * (n_material - 1))

    # Save sag map for Zemax import
    np.savetxt('metalens_sag_LCP_for_zemax.txt', sag_map, fmt='%.8e')

    print("Sag map for LCP saved as 'metalens_sag_LCP_for_zemax.txt'")

# Optional: Plot the unwrapped phase
plt.figure()
plt.pcolormesh(x * 1e6, y * 1e6, unwrapped_phase, shading='auto', cmap='jet')
plt.colorbar(label='Phase (rad)')
plt.xlabel('x (μm)')
plt.ylabel('y (μm)')
plt.title('Unwrapped Phase at Metalens Exit (LCP)')
plt.axis('equal')
plt.show()

# Optional: Plot the sag map
plt.figure()
plt.pcolormesh(x * 1e6, y * 1e6, sag_map, shading='auto', cmap='jet')
plt.colorbar(label='Sag (m)')
plt.xlabel('x (μm)')
plt.ylabel('y (μm)')
plt.title('Metalens Surface Sag for Zemax (LCP)')
plt.axis('equal')
plt.show()