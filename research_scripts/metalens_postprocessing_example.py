import numpy as np
# from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import importlib.util
import sys
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
# lumfilename = "metalens_from_p1000_d_f_T_phase_R_A_LCP.txt_H_on_RCP_R15.0um_f10.0um_Marco approach_2025-02-21_16-33-54.fsp"

magnetic_field_ON = 1 # 0 for H off # 1 for H on

circ_pol = 1 # 0 for linear # 1 for circular
##here should change
LCP = 0
RCP = 1-LCP

LUT_from_COMSOL = 1
LUT_from_LUMERICAL = 1-LUT_from_COMSOL

use_LUT_no_H = 0 # 1 to use LUT calculated for no H, 0 to use LUT corresponding to the excitation configuration
use_pump_as_in_LUT = 0 # 1 to use H and LCP or RCP as in comsol filename, 0 to use ones stated above

make_metalens_Marco_approach = 1 # 1 to use Marco's approach, 0 - to use Lumerical's approach
make_metalens_Lumerical_approach = 1 - make_metalens_Marco_approach  # Ensures only one is 1

wavelength = 1.55e-6
period = 1000e-9
height = 1300e-9
mat_pillar = np.bytes_("Object defined dielectric")
mat_sub = np.bytes_("Object defined dielectric")
index_pillar = 1 # material BIG is introduced later
index_sub = 1.4491
n_subs = 1.4491
E0 = 1
f = 10e-6  # Focal length in meters
R = 15e-6  # Lens radius in meters
monitor_name2 = "DFTmonitorXZ"
monitor_name = "DFTmonitor0"
lens_radius = R

if make_metalens_Lumerical_approach==1:
    plot_title123 = 'Lumerical approach'
elif make_metalens_Marco_approach==1:
    plot_title123 = 'Marco approach'

if LUT_from_COMSOL ==1:
    # Load the data from the text file, skipping comment lines (starting with '%')
    if LCP == 1:
        LookUp_table_filename = 'p1000_d_f_T_phase_R_A_LCP_NO_H_fitto.txt'
        your_input_file_is_from_ = -1 #
        if magnetic_field_ON == 1:
            LookUp_table_filename = 'p1000_d_f_T_phase_R_A_LCP.txt'
            your_input_file_is_from_ = -2 #
    if RCP == 1:
        # plot_title3 = 'RCP'
        LookUp_table_filename = 'p1000_d_f_T_phase_R_A_RCP.txt'
        your_input_file_is_from_ = -3 #

    if use_LUT_no_H == 1:
        LookUp_table_filename = 'p1000_d_f_T_phase_R_A_LCP_NO_H_fitto.txt'
        your_input_file_is_from_ = -1 #

    if use_pump_as_in_LUT ==1:
        if your_input_file_is_from_ == -1:
            circ_pol = 1 # 0 for linear # 1 for circular
            LCP = 1
            RCP = 0
            magnetic_field_ON = 0
        if your_input_file_is_from_ == -2:
            circ_pol = 1 # 0 for linear # 1 for circular
            LCP = 1
            RCP = 0
            magnetic_field_ON = 1
        if your_input_file_is_from_ == -3:
            circ_pol = 1 # 0 for linear # 1 for circular
            LCP = 0
            RCP = 1
            magnetic_field_ON = 1
    else: # define pump polarization here
        LCP = 0
        RCP = 1-LCP




    if your_input_file_is_from_ < 0: # Check the input file
        # If the file has comments ('%' in MATLAB), use a more flexible approach to ignore them
        data_list = []
        with open(LookUp_table_filename, 'r') as file:
            for line in file:
                if not line.startswith('%'):  # Ignore comment lines
                    data_list.append(line.strip().split())  # Split each line into columns

    # Convert list to NumPy array
    dataComsol = np.array(data_list, dtype=float)
    radius = dataComsol[:, 0]/2
    phase_comsol = (dataComsol[:, 3])
    r_min = np.min(radius)
    r_max = np.max(radius)
    # Unwrap and normalize the phase data
    phase_comsol_unwrap = np.unwrap(phase_comsol)
    phase_comsol0 = phase_comsol[0]
    phase = np.abs(phase_comsol_unwrap) - np.abs(phase_comsol0)

np_points = round(R / period)
x_mask = np.arange(-np_points, np_points + 1) * period
y_mask = x_mask
phase_shift = 0*np.pi# 0 # 1.5
iph = np.argmin(np.abs(phase - phase_shift))
phase_target = np.mod(((-(2 * np.pi / wavelength) * (np.sqrt(x_mask**2 + f**2) - f) ) + phase[iph]), 2 * np.pi)
radius_mask = np.array([radius[np.argmin(np.abs(phase - pt))] for pt in phase_target])


if magnetic_field_ON == 1:
    plot_title1 = 'Magnetic_field_ON'
else: plot_title1 = 'Magnetic_field_off'

if circ_pol == 1:
    plot_title2 = 'Circular_polarization'
    if LCP == 1:
        plot_title3 = 'LCP'
    if RCP == 1:
        plot_title3 = 'RCP'
    plot_title11 = f'{plot_title2}_{plot_title3}'
else: plot_title11 = 'Linear_polarization'

plot_title = f'{plot_title1}_{plot_title11}_FDTD_{plot_title123}_1'


# Open Lumerical FDTD
with lumapi.FDTD() as fdtd:
    fdtd.load(lumfilename)
    monitor_name2 = "DFTmonitorXZ"
    print(fdtd.getresult("DFTmonitorXZ"))  # Check available data fields
    # # # # # # # # # # # # # # # # #
    # Postprocessing
    # # # # # # # # # # # # # # # # #
    index_sub = n_subs
    # # Get refractive index of the substrate and of the air
    # n2 = float(fdtd.getnamed("FDTD", "index"))
    # # n1_pil = float(fdtd.getnamed("metalens::circle", "index"))
    # n1 = float(fdtd.getnamed("substrate", "index"))
    n1 = index_sub # index_sub
    n2 = index_sub
    rho = (n1 - n2) / (n1 + n2)
    Tsa = - abs(rho)**2  # Transmittance from the substrate to the air, without the metalens

    # Intensity of the incident light
    # E0 = 1  # [V/m]
    Z0_const = 376.73
    I0 = 1 / (2 * Z0_const) * n1 * E0**2  # Intensity of the incident light
    Inc_Pow_aperture = I0 * R**2 * np.pi  # Power entering in the circular aperture

    # Near-field

    print(fdtd.getresult("DFTmonitorXZ"))  # Check available data fields

    ExzGB = fdtd.getresult("DFTmonitorXZ", "E")
    xGB = ExzGB['x']
    xGB = xGB[:, 0]
    zGB = ExzGB['z']
    zGB = zGB[:, 0]
    # # # # #        coord   wwl components(Ex, Ey,...)
    # #               x  y  z  wwl x=0 y=1 z=2
    ExGB = ExzGB['E'][:, 0, :, 0, 0]
    EyGB = ExzGB['E'][:, 0, :, 0, 1]
    EzGB = ExzGB['E'][:, 0, :, 0, 2]

    EGB2 = (np.abs(ExGB + EyGB + EzGB))**2

    # Create a 1x3 grid of subplots (1 row, 3 columns)
    # fig, axes = plt.subplots(1, 3, figsize=(8, 3))
    # fig.suptitle(plot_title)
    # cf1 = axes[0].contourf(zGB* 1e6, xGB* 1e6, np.real(ExGB), 100)
    # axes[0].set_title('Ex')
    # axes[0].set_xlabel('x (um)')
    # axes[0].set_ylabel('z (um)')
    # axes[0].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf1, ax=axes[0], label="Enf [V/m]")
    # cf2 = axes[1].contourf(zGB* 1e6, xGB* 1e6, np.real(EyGB), 100)
    # axes[1].set_title('Ey')
    # axes[1].set_xlabel('x (um)')
    # axes[1].set_ylabel('z (um)')
    # axes[1].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf2, ax=axes[1], label="Enf [V/m]")
    # cf3 = axes[2].contourf(zGB* 1e6, xGB* 1e6, np.real(EzGB), 100)
    # axes[2].set_title('Ez')
    # axes[2].set_xlabel('x (um)')
    # axes[2].set_ylabel('z (um)')
    # axes[2].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf3, ax=axes[2], label="Enf [V/m]")
    # # Adjust layout to prevent overlap
    # plt.tight_layout()
    # # Show the plot
    # plt.show()
    # plt.contourf(zGB* 1e6, xGB* 1e6, ExGB, 100)
    # plt.colorbar()
    # plt.title('Ex')
    # plt.ylabel('x (um)')
    # plt.xlabel('z (um)')
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()



    # Create a 1x3 grid of subplots (1 row, 3 columns)
    # fig, axes = plt.subplots(1, 3, figsize=(8, 3))
    # fig.suptitle(plot_title)
    # cf1 = axes[0].contourf(zGB* 1e6, xGB* 1e6, np.angle((ExGB)), 100)
    # axes[0].set_title('phase Ex')
    # axes[0].set_xlabel('x (um)')
    # axes[0].set_ylabel('z (um)')
    # axes[0].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf1, ax=axes[0], label="phase")
    # cf2 = axes[1].contourf(zGB* 1e6, xGB* 1e6, np.angle((EyGB)), 100)
    # axes[1].set_title('phase Ey')
    # axes[1].set_xlabel('x (um)')
    # axes[1].set_ylabel('z (um)')
    # axes[1].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf2, ax=axes[1], label="phase")
    # cf3 = axes[2].contourf(zGB* 1e6, xGB* 1e6, np.angle((EzGB)), 100)
    # axes[2].set_title('phase Ez')
    # axes[2].set_xlabel('x (um)')
    # axes[2].set_ylabel('z (um)')
    # axes[2].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf3, ax=axes[2], label="phase")
    # # Adjust layout to prevent overlap
    # plt.tight_layout()
    # # Show the plot
    # plt.show()

    # plt.contourf(zGB* 1e6, xGB* 1e6,  EGB2, 100) # /np.pi
    # plt.colorbar()
    # plt.title('|Enf|^2')
    # plt.xlabel('x (um)')
    # plt.ylabel('z (um)')
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()

    EyzGB1 = fdtd.getresult("DFTmonitorYZ", "E")
    yGB1 = EyzGB1['y']
    yGB1 = yGB1[:, 0]
    zGB1 = EyzGB1['z']
    zGB1 = zGB1[:, 0]
    # # # # #        coord   wwl components(Ex, Ey,...)
    # #               x  y  z  wwl x=0 y=1 z=2
    ExGB1 = EyzGB1['E'][0, :, :, 0, 0]
    EyGB1 = EyzGB1['E'][0, :, :, 0, 1]
    EzGB1 = EyzGB1['E'][0, :, :, 0, 2]

    EGB1 = (np.abs(ExGB1 + EyGB1 + EzGB1))**2


    # plt.contourf(zGB1* 1e6, yGB1* 1e6,  EGB1, 100) # /np.pi
    # plt.colorbar()
    # plt.title('|Enf|^2')
    # plt.xlabel('y (um)')
    # plt.ylabel('z (um)')
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    Enf = fdtd.getresult("DFTmonitor0", "E")
    Hnf = fdtd.getresult("DFTmonitor0", "H")
    Lambda = Enf['lambda']
    Lambda = Lambda[:, 0]


    x = Enf['x']
    x = x[:, 0]
    y = Enf['y']
    y = y[:, 0]
    # # # # #        coord   wwl components(Ex, Ey,...)
    # #             x  y  z  wwl x=0 y=1 z=2
    Enfx = Enf['E'][:, :, 0, :, 0]
    Enfy = Enf['E'][:, :, 0, :, 1]
    Enfz = Enf['E'][:, :, 0, :, 2]
    Hnfx = Hnf['H'][:, :, 0, :, 0]
    Hnfy = Hnf['H'][:, :, 0, :, 1]

    Enf_xy = Enfx[:, :, 0]+Enfy[:, :, 0]+Enfz[:, :, 0]

    Pznf = np.real(Enfx * np.conj(Hnfy) - Enfy * np.conj(Hnfx)) / 2
    # Power_on_monitor_nf = integrate(Pznf, slice(1, 2), x, y)
    Power_on_monitor_nf = np.trapz(np.trapz(Pznf, x, axis=1), y, axis=0)

    Tml = Power_on_monitor_nf / Inc_Pow_aperture  # Power arriving on the monitor field divided by the power of the plane wave passing through the slit


    # Create a 1x3 grid of subplots (1 row, 3 columns)
    # fig, axes = plt.subplots(1, 3, figsize=(10, 3))
    # fig.suptitle(plot_title)
    # cf1 = axes[0].contourf(x* 1e6, y* 1e6, np.real(Enfx[:, :, 0]), 100)
    # axes[0].set_title('Ex')
    # axes[0].set_xlabel('x (um)')
    # axes[0].set_ylabel('y (um)')
    # axes[0].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf1, ax=axes[0], label="Enf [V/m]")
    # cf2 = axes[1].contourf(x* 1e6, y* 1e6, np.real(Enfy[:, :, 0]), 100)
    # axes[1].set_title('Ey')
    # axes[1].set_xlabel('x (um)')
    # axes[1].set_ylabel('y (um)')
    # axes[1].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf2, ax=axes[1], label="Enf [V/m]")
    # cf3 = axes[2].contourf(x* 1e6, y* 1e6, np.real(Enfz[:, :, 0]), 100)
    # axes[2].set_title('Ez')
    # axes[2].set_xlabel('x (um)')
    # axes[2].set_ylabel('y (um)')
    # axes[2].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf3, ax=axes[2], label="Enf [V/m]")
    # # Adjust layout to prevent overlap
    # plt.tight_layout()
    # # Show the plot
    # plt.show()
    #
    # Create a 1x3 grid of subplots (1 row, 3 columns)
    # fig, axes = plt.subplots(1, 3, figsize=(10, 3))
    # fig.suptitle(LookUp_table_filename)
    # cf1 = axes[0].contourf(x* 1e6, y* 1e6, np.angle(Enfx[:, :, 0]), 100) # , vmin=-np.pi, vmax=np.pi)
    # axes[0].set_title('phase Ex')
    # axes[0].set_xlabel('x (um)')
    # axes[0].set_ylabel('y (um)')
    # axes[0].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf1, ax=axes[0], label="phase")
    # cf2 = axes[1].contourf(x* 1e6, y* 1e6, np.angle(Enfy[:, :, 0]), 100)# , vmin=-np.pi, vmax=np.pi)
    # axes[1].set_title('phase Ey')
    # axes[1].set_xlabel('x (um)')
    # axes[1].set_ylabel('y (um)')
    # axes[1].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf2, ax=axes[1], label="phase")
    # cf3 = axes[2].contourf(x* 1e6, y* 1e6, np.angle(Enfz[:, :, 0]), 100)# , vmin=-np.pi, vmax=np.pi)
    # axes[2].set_title('phase Ez')
    # axes[2].set_xlabel('x (um)')
    # axes[2].set_ylabel('y (um)')
    # axes[2].set_aspect('equal', adjustable='box')
    # fig.colorbar(cf3, ax=axes[2], label="phase")
    # # Adjust layout to prevent overlap
    # plt.tight_layout()
    # # Show the plot
    # plt.show()


    # plt.contourf(x* 1e6, y* 1e6, np.angle((Enf_xy)), 100) # /np.pi
    # plt.colorbar()
    # plt.title('phase(E)')
    # plt.xlabel('x (um)')
    # plt.ylabel('y (um)')
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    #
    # Phase of the near-field
    phase_x = np.mod(np.unwrap(np.angle((Enfx[:, (len(y) // 2), 0]))), 2*np.pi)
    phase_y = np.mod(np.unwrap(np.angle((Enfy[:, (len(y) // 2), 0]))), 2*np.pi)
    phase_nf = np.mod(np.unwrap(np.angle(Enf_xy[:, (len(y) // 2)])), 2*np.pi)

    # phase_x = phase_x - (np.max(phase_x) - 2 * np.pi)
    phase_target1 = (2 * np.pi / wavelength) * (f - np.sqrt(x**2 + f**2))  # Hyperbolic lens is assumed

    # Plot the phase
    # plt.figure()
    # plt.bar(x_mask * 1e6, phase_target, width=period * 1e6, color='red', alpha=0.7, label="Target",edgecolor='black')
    # # plt.step(x_mask * 1e6, phase_target, where='mid', color='black', marker='.')
    # plt.plot(x_mask * 1e6, phase_target, label="Target Phase vs. Position", color='blue', marker='o')
    # plt.plot(x * 1e6, phase_x, label="Simulation", marker='.')
    # plt.xlabel("Position (μm)")
    # plt.ylabel("Phase Ex (rad)")
    # plt.title("Phase Ex (rad)")
    # plt.legend()
    # plt.xlim([0, (R + period / 2) * 1e6])
    # plt.grid(True)
    # plt.show()
    #
    # plt.plot(x * 1e6, np.mod(phase_target1, 2*np.pi), label="Target Phase approx", marker='.') # linestyle='None',
    # plt.plot(x_mask * 1e6, np.mod(phase_target, 2*np.pi), label="Target", color='red', marker='*')
    # plt.plot(x * 1e6, phase_x, label="Simulation", marker='.')
    # plt.xlabel("x (um)")
    # plt.ylabel("Phase Ex (rad)")
    # plt.legend()
    # plt.xlim([-R * 1e6, R * 1e6])
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    #
    # plt.plot(x_mask * 1e6, np.mod(phase_target, 2*np.pi), label="Target", color='red', marker='*')
    # plt.plot(x * 1e6, phase_y, label="Simulation", marker='.')
    # plt.xlabel("x (um)")
    # plt.ylabel("Phase Ey (rad)")
    # plt.legend()
    # plt.xlim([-R * 1e6, R * 1e6])
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    # Plot the phase
    # plt.figure()
    # plt.bar(x_mask * 1e6, phase_target, width=period * 1e6, color='red', alpha=0.7, label="Target",edgecolor='black')
    # plt.plot(x * 1e6, phase_y, label="Simulation", marker='.')
    # plt.xlabel("Position (μm)")
    # plt.ylabel("Phase Ey (rad)")
    # plt.title("Phase Ey (rad)")
    # plt.legend()
    # plt.xlim([0, (R + period / 2) * 1e6])
    # plt.grid(True)
    # plt.show()
    # # Plot the phase
    # plt.plot(x * 1e6, np.mod(phase_target1, 2*np.pi), label="Target Phase approx", marker='.') # linestyle='None',
    # plt.plot(x_mask * 1e6, np.mod(phase_target, 2*np.pi), label="Target", color='red', marker='*')
    # plt.plot(x * 1e6, phase_nf, label="Simulation", marker='.')
    # plt.xlabel("x (um)")
    # plt.ylabel("Phase E (rad)")
    # plt.legend()
    # plt.xlim([-R * 1e6, R * 1e6])
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    #
    # Plot the phase
    # plt.figure()
    # plt.bar(x_mask * 1e6, phase_target, width=period * 1e6, color='red', alpha=0.7, label="Target",edgecolor='black')
    # plt.plot(x * 1e6, phase_nf, label="Simulation", marker='.')
    # plt.xlabel("Position (μm)")
    # plt.ylabel("Phase E (rad)")
    # plt.title("Phase E (rad)")
    # plt.legend()
    # plt.xlim([0, (R + period / 2) * 1e6])
    # plt.grid(True)
    # plt.show()
    #

    Enf1 = fdtd.getresult("DFTmonitor1", "E")
    x1 = Enf1['x']
    x1 = x1[:, 0]
    y1 = Enf1['y']
    y1 = y1[:, 0]
    # # # # #        coord   wwl components(Ex, Ey,...)
    # #             x  y  z  wwl x=0 y=1 z=2
    Enfx1 = Enf1['E'][:, :, 0, :, 0]
    Enfy1 = Enf1['E'][:, :, 0, :, 1]
    Enfz1 = Enf1['E'][:, :, 0, :, 2]
    Enf_xy1 = Enfx1[:, :, 0]+Enfy1[:, :, 0]+Enfz1[:, :, 0]
    zmonitor1 = f + height/2
    # plt.contourf(x1 * 1e6, y1 * 1e6, np.abs((Enf_xy1))**2, 100) # /np.pi
    # plt.colorbar()
    # plt.title(f'|Enf|^2_z=_{zmonitor1*1e6}_um')
    # plt.xlabel('x (um)')
    # plt.ylabel('y (um)')
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    #

    # # ######################
    # # # Farfield
    # # ######################

    # # Downsample the near field to speed up far-field projections
    fdtd.farfieldsettings("override near field mesh", True)
    fdtd.farfieldsettings("near field samples per wavelength", 2)  # At Nyquist limit

    # Choose area to plot and spatial resolution
    # x_val = np.linspace(-12e-6, 12e-6, 400)
    # y_val = np.linspace(-12e-6, 12e-6, 400)
    # z_val = np.linspace(-5e-6, 12e-6, 400)
    x_val = np.linspace(-1.1*R, 1.1*R, 401)
    y_val = np.linspace(-1.1*R, 1.1*R, 401)
    z_val = np.linspace(-2*height, 3*f, 401)

    # Plot E field intensity profile as a function of xz and xy (focal) planes
    E_xz = fdtd.farfieldexact3d(monitor_name, x_val, 0, z_val)
    E2_xz = np.abs(E_xz[:, 0, :, 0])**2 + np.abs(E_xz[:, 0, :, 1])**2 + np.abs(E_xz[:, 0, :, 2])**2

    # plt.imshow(E2_xz, aspect="auto", extent=[x_val.min() * 1e6, x_val.max() * 1e6, z_val.min() * 1e6, z_val.max() * 1e6])
    # plt.xlabel("z (um)")
    # plt.ylabel("x (um)")
    # # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.title("|E|^2 (x-z plane)")
    # plt.colorbar()
    # plt.show()

    # plt.contourf(z_val* 1e6, x_val* 1e6, np.real((E2_xz)), 100)
    # plt.colorbar()
    # plt.title('Eff^2')
    # plt.xlabel('z (um)')
    # plt.ylabel('x (um)')
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    #
    # Far-field along z-axis
    index_x0 = np.argmin(np.abs(x_val))#  np.where(x_val == 0)[0]
    E2_z = E2_xz[index_x0, :]

    # Plot far-field along z-axis
    # plt.plot(z_val * 1e6, E2_z)
    # plt.xlabel("z (um)")
    # plt.ylabel("|Eff|^2 (V^2/m^2)")
    # plt.title("Farfield along z-axis")
    # plt.show()
    #
    # # Far-field along z-axis
    # E_z = fdtd.farfieldexact3d(monitor_name, 0, 0, z_val)
    # E2_z = np.abs(E_z[0, 0, :, 0])**2 + np.abs(E_z[0, 0, :, 1])**2 + np.abs(E_z[0, 0, :, 2])**2 # Sum along the 4th index (Ex^2 + Ey^2 + Ez^2)

    # # Plot far-field along z-axis
    # plt.plot(z_val * 1e6, E2_z)
    # plt.xlabel("z (um)")
    # plt.ylabel("|E|^2 (V^2/m^2)")
    # plt.title("Farfield along z-axis")
    # plt.show()

    # Find focal plane
    ind = np.argmax(E2_z)
    zm = z_val[ind]

    # Far-field in x-direction at the focal point (x=0, y=0, z=zm)
    E_x = fdtd.farfieldexact3d(monitor_name, x_val, 0, zm)
    E2_x = np.abs(E_x[:, 0, 0, 0])**2 + np.abs(E_x[:, 0, 0, 1])**2 + np.abs(E_x[:, 0, 0, 2])**2 # Sum along the 4th index (Ex^2 + Ey^2 + Ez^2)

    # plt.plot(x_val * 1e6, E2_x)
    # plt.xlabel("x (um)")
    # plt.ylabel("|Eff|^2 (V^2/m^2)")
    # plt.title(f'Farfield at focal plane z={zm*1e6:.2f} um')
    # plt.show()
    #
    # Far-field in y-direction at the focal point (x=0, y=0, z=zm)
    E_y = fdtd.farfieldexact3d(monitor_name, 0, y_val, zm)
    E2_y = np.abs(E_y[0, :, 0, 0])**2 + np.abs(E_y[0, :, 0, 1])**2 + np.abs(E_y[0, :, 0, 2])**2 # Sum along the 4th index (Ex^2 + Ey^2 + Ez^2)

    # plt.plot(y_val * 1e6, E2_y)
    # plt.xlabel("y (um)")
    # plt.ylabel("|Eff|^2 (V^2/m^2)")
    # # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.title("Farfield at focal plane (y-direction)")
    # plt.show()
    #


    E_xy = fdtd.farfieldexact3d(monitor_name, x_val, y_val, zm)
    E__xy = (E_xy[:, :, 0, 0]) + (E_xy[:, :, 0, 1]) + (E_xy[:, :, 0, 2])
    E2_xy = np.abs(E_xy[:, :, 0, 0])**2 + np.abs(E_xy[:, :, 0, 1])**2 + np.abs(E_xy[:, :, 0, 2])**2

    # plt.imshow(E2_xy, aspect="auto", extent=[x_val.min() * 1e6, x_val.max() * 1e6, y_val.min() * 1e6, y_val.max() * 1e6])
    # plt.xlabel("x (um)")
    # plt.ylabel("y (um)")
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.title("|Eff|^2 (x-y plane)")
    # plt.colorbar()
    # plt.show()
    #
    # # ##################
    # # # Part 2: Focusing efficiency with Poynting vector
    # # ##################

    f_monitor = fdtd.getdata(monitor_name, "f")

    x_vals = np.linspace(-20 * 1e-6, 20 * 1e-6, 41)
    y_vals = np.linspace(-20 * 1e-6, 20 * 1e-6, 41)
    X, Y = np.meshgrid(x_vals, y_vals)

    EHff = fdtd.farfieldexact3d("DFTmonitor0", x_vals, y_vals, zm, {"f": 1, "field": "E and H", "index": 1})


    Exff = EHff['E'][:, :, 0, :, 0]
    Eyff = EHff['E'][:, :, 0, :, 1]
    Ezff = EHff['E'][:, :, 0, :, 2]
    Hxff = EHff['H'][:, :, 0, :, 0]
    Hyff = EHff['H'][:, :, 0, :, 1]
    Hzff = EHff['H'][:, :, 0, :, 2]


    # # Calculate components of Poynting vector
    Pzff = np.real(Exff * np.conj(Hyff) - Eyff * np.conj(Hxff)) / 2

    # Filter: Circle with R = second Airy zero
    secondAiryZero = 1.22 * (2 * wavelength * f) / ((2*R)) # x_vals - 7.0156 / ((2 * np.pi / wavelength) * (R / f))
    ind1 = np.argmin(abs(x_vals - secondAiryZero))
    # ind1 = np.where(abs(secondAiryZero) == 0)
    # ind1 = np.where(abs(secondAiryZero) <7e-8)
    filter_ = (X**2 + Y**2) <= x_vals[ind1]**2

    # Integrate Pz over the filter region
    Pzff2d = np.squeeze(Pzff)  # Remove unnecessary dimensions
    Power_around_focus = np.trapz(np.trapz(np.multiply(Pzff2d, filter_), x_vals, axis=0), y_vals, axis=0)

    efficiency_on_the_transmitted = Power_around_focus / Power_on_monitor_nf
    efficiency_on_the_incident = Power_around_focus / Inc_Pow_aperture

    # Plot filter and transmission
    # plt.imshow(filter_, extent=[x_vals.min() * 1e6, x_vals.max() * 1e6, y_vals.min() * 1e6, y_vals.max() * 1e6])
    # plt.xlabel("x (um)")
    # plt.ylabel("y (um)")
    # plt.title("Integration Region")
    # plt.colorbar()
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    #
    # c = 3e8  # Speed of light in m/s
    # plt.plot(wavelength * 1e6, efficiency_on_the_transmitted * 100, label="Transmitted Efficiency", linestyle='None', marker='*')
    # plt.plot(wavelength * 1e6, efficiency_on_the_incident * 100, label="Incident Efficiency", linestyle='None', marker='*')
    # plt.xlabel("lambda (um)")
    # plt.ylabel("Focusing Efficiency (%)")
    # # plt.ylim(0, 100)
    # plt.grid(True)
    # plt.legend()
    # plt.show()
    #
    # Find the maximum intensity and half-max value
    max_intensity = np.max(E2_xy)
    half_max = max_intensity / 2
    # Pick a row (e.g., the center row) for the x-direction FWHM
    center_y_index = E2_xy.shape[0] // 2  # Middle row
    intensity_profile_x = E2_xy[center_y_index, :]  # Intensity profile along x-axis
    # Find the indices where the intensity exceeds half-max in the x-direction
    fwhm_indices_x = np.where(intensity_profile_x >= half_max)[0]
    # Calculate the FWHM in the x-direction
    FWHM_x = (x_val[fwhm_indices_x[-1]] - x_val[fwhm_indices_x[0]]) * 1e6  # Convert to micrometers
    # Pick a column (e.g., the center column) for the y-direction FWHM
    center_x_index = E2_xy.shape[1] // 2  # Middle column
    intensity_profile_y = E2_xy[:, center_x_index]  # Intensity profile along y-axis
    # Find the indices where the intensity exceeds half-max in the y-direction
    fwhm_indices_y = np.where(intensity_profile_y >= half_max)[0]
    # Calculate the FWHM in the y-direction
    FWHM_y = (y_val[fwhm_indices_y[-1]] - y_val[fwhm_indices_y[0]]) * 1e6  # Convert to micrometers
    print(f"fwhm_indices_x: {fwhm_indices_x} ")
    print(f"fwhm_indices_y: {fwhm_indices_y} ")
    print(f"FWHM in x-direction: {FWHM_x} µm")
    print(f"FWHM in y-direction: {FWHM_y} µm")


    phase_at_focus = np.angle(E__xy)  # Phase of the wavefront
    # Unwrap phase and subtract the mean to center it
    aberration_map = np.unwrap(phase_at_focus) - np.mean(np.unwrap(phase_at_focus))
    # Visualize the aberration map
    # plt.contourf(x_val* 1e6, y_val* 1e6, aberration_map, 100)
    # plt.colorbar()
    # plt.title('Wavefront Aberration Map')
    # plt.xlabel('x (um)')
    # plt.ylabel('y (um)')
    # plt.axis('equal')  # Ensure aspect ratio is equal
    # plt.show()
    #
    # # from pyzernike import Zernike
    # from pyzernike import polynomials
    from scipy.special import factorial
    from skimage.restoration import unwrap_phase
    import math

    # import zernike
    # print(dir(zernike))


    # import pyzernike
    # print(dir(pyzernike))

    def zernike_radial(n, m, r, raggio):
        """Calculate the radial Zernike polynomial R_n^m(r)."""
        R = np.zeros_like(r, dtype=float)
        for k in range((n - m) // 2 + 1):
            num = (-1)**k * factorial(n - k)
            den = factorial(k) * factorial((n + m) / 2 - k) * factorial((n - m) / 2 - k)
            R += num / den * (r/raggio)**(n - 2 * k)
        if m!=0:
            R=R/math.sqrt(math.pi/(2*n+2))
        else:
            R=R/math.sqrt(2*math.pi/(2*n+2))
        return R

    def zernike_polynomial(n, m, r, theta, raggio):
        """Calculate the full Zernike polynomial Z_n^m(r, theta)."""
        R = zernike_radial(n, abs(m), r, raggio)
        if m >= 0:
            return R * np.cos(m * theta)
        elif m < 0:
            return R * np.sin(-m * theta)

    def fit_zernike(wavefront, r, theta, n_max,dx,dy,raggio):
        """Fit Zernike polynomials to a wavefront and return coefficients."""
        coefficients = []
        matrix_coefficients=np.empty((n_max+1,2*n_max+1,))
        matrix_coefficients[:] = np.nan
        for n in range(n_max + 1):
            for m in range(-n, n + 1, 2):  # m steps by 2 to satisfy Zernike conditions
                Z_nm = zernike_polynomial(n, m, r, theta, raggio)
                c_nm = np.sum(wavefront*Z_nm)*(dx)*(dy)/raggio**2# Weighted least-squares
                coefficients.append((n, m, c_nm))
                matrix_coefficients[n][n_max+m]=c_nm

        n_lst= np.arange(0, n_max+1)
        m_lst= np.arange(-n_max, n_max+1)
        return coefficients,matrix_coefficients,n_lst,m_lst

    radius_of_the_lens = R
    lda = wavelength
    focal_of_the_lens = f

    # x_grid = np.linspace(-radius_of_the_lens, radius_of_the_lens, 1000)
    # y_grid = np.linspace(-radius_of_the_lens, radius_of_the_lens, 1000)

    # x_size = len(x_grid)
    # y_size = len(y_grid)
    # E_x_array=np.zeros((x_size, y_size), dtype=complex)

    # # # z=2.3e-6 # Why did you choose this coordinate?
    # # # # Find the index of the closest point to (x, y, z)
    # # # z_index = np.argmin(np.abs(z_lum - z))
    # z_index = ind
    # for iii in range(x_size):

    #     x_index = np.argmin(np.abs(x_val - x_grid[iii]))

    #     for jjj in range(y_size):

    #         y_index = np.argmin(np.abs(y_val - y_grid[jjj]))
    #         E_x_array[iii][jjj]= data[x_index, y_index, z_index, :][0]
    x_grid = x_val
    y_grid = y_val
    E_x_array = E__xy #now it's actually Ex+Ey+Ez
    X, Y = np.meshgrid(x_grid, y_grid)


    RR=np.sqrt(X**2+Y**2)
    Theta = np.arctan2(Y, X)

    # Plotting
    # plt.figure(figsize=(8, 6))
    # plt.pcolormesh(X*1e6, Y*1e6, np.abs(E_x_array), shading='auto', cmap='jet')  # Create a 2D color plot
    # plt.colorbar(label='|E|')
    # # plt.colorbar(label='Field Intensity')  # Add a colorbar
    # plt.xlabel('x (um)')
    # plt.ylabel('y (um)')
    # plt.title('2D Matrix Plot')
    # plt.show()



    #Extract phase
    raw_phase = np.angle(E_x_array)  # Modulo 2pi phase
    unwrapped_phase = unwrap_phase(raw_phase)  # Unwrap the phase

    nX=X.shape
    for iii in range(nX[0]-1):
        for jjj in range(nX[1]-1):
            if RR[iii][jjj]/radius_of_the_lens>1:
                unwrapped_phase[iii][jjj]=0




    wavefront = lda*unwrapped_phase/2/math.pi;
    wavefront_p=-RR**2/2/focal_of_the_lens

    for iii in range(nX[0]-1):
        for jjj in range(nX[1]-1):
            if RR[iii][jjj]/radius_of_the_lens>1:
                wavefront_p[iii][jjj]=0

    dW=wavefront-wavefront_p



    for iii in range(nX[0]-1):
        for jjj in range(nX[1]-1):
            if RR[iii][jjj]/radius_of_the_lens>1:
                dW[iii][jjj]=0
    # Plotting
    # plt.figure
    # plt.pcolormesh(X*1e6, Y*1e6, raw_phase, shading='auto', cmap='jet')  # Create a 2D color plot
    # plt.colorbar(label='unwr. phase')  # Add a colorbar
    # plt.xlabel('x (um)')
    # plt.ylabel('y (um)')
    # plt.axis('equal')
    # plt.title('raw phase')
    # plt.show()
    #
    # Plotting
    # plt.figure(figsize=(8, 6))
    # plt.pcolormesh(X*1e6, Y*1e6, unwrapped_phase, shading='auto', cmap='jet')  # Create a 2D color plot
    # plt.colorbar(label='unwr. phase')  # Add a colorbar
    # plt.xlabel('x (um)')
    # plt.ylabel('y (um)')
    # plt.title('unwrapped phase')
    # plt.show()

    ##################
    # plt.figure(figsize=(8, 6))
    # plt.pcolormesh(X*1e6, Y*1e6, dW, shading='auto', cmap='jet')  # Create a 2D color plot
    # plt.colorbar(label='')  # Add a colorbar
    # plt.xlabel('m')
    # plt.ylabel('n')
    # plt.title('dW')
    # plt.show()
    # plt.rcParams.update({'font.size': 22})


    # Fit Zernike coefficients
    # n_max = 4  # Maximum order
    # [zernike_coeffs,matr_coeff,n_list,m_list] = fit_zernike(dW, RR, Theta, n_max,x_grid[1]-x_grid[0],y_grid[1]-y_grid[0],radius_of_the_lens)
    #
    # # Print Zernike coefficients
    # print("Zernike Coefficients (n, m, c_nm):")
    # for coeff in zernike_coeffs:
    #     print(coeff)



    #plot Zernike Coefficients
    #ML,NL = np.meshgrid(m_list,n_list)

    #plt.figure(figsize=(8, 6))
    #plt.pcolormesh(ML, NL, np.abs(matr_coeff), shading='auto', cmap='jet')  # Create a 2D color plot
    #plt.gca().invert_yaxis()
    #plt.colorbar(label='weight (a.u.)')  # Add a colorbar
    #plt.xlabel('m')
    #plt.ylabel('n')
    #plt.title('Zernike pol.')
    #plt.show()
    #plt.rcParams.update({'font.size': 22})



    #
    # Save to .mat file
    # # # # # # # # # # #
    import scipy.io
    # lens_design = f'{lumfilename}_1.mat'
    # scipy.io.savemat(lens_design, {
    #     "position_m": x_mask,
    #     "phase_target": phase_target,
    #     "radius_pillar": radius_mask,
    #     "x_ff": x_val,
    #     "y_ff": y_val,
    #     "z_ff": z_val,
    #     "E2_x_ff": E2_x,
    #     "E2_y_ff": E2_y,
    #     "E2_z_ff": E2_z,
    # }, do_compression=True)
    # print(f"LUT data saved to {lens_design}")

    lens_result = f'metalens_result_{lumfilename}.mat'
    scipy.io.savemat(lens_result, {
        "position_m": x_mask,
        "phase_target": phase_target,
        "radius_pillar": radius_mask,
        "wavelength": wavelength,
        "height": height,
        "period": period,
        "mat_pillar": mat_pillar,
        "index_pillar": index_pillar,
        "mat_sub": mat_sub,
        "index_sub": index_sub,
        "phase": phase,
        "radius": radius,
        "x_xz_nf": xGB,
        "z_xz_nf": zGB,
        "Ex_xz_nf": ExGB,
        "Ey_xz_nf": EyGB,
        "Ez_xz_nf": EzGB,
        "y_yz_nf": yGB1,
        "z_yz_nf": zGB1,
        "Ex_yz_nf": ExGB1,
        "Ey_yz_nf": EyGB1,
        "Ez_yz_nf": EzGB1,
        "x_xy_nf": x,
        "y_xy_nf": y,
        "Ex_xy_nf": Enfx,
        "Ey_xy_nf": Enfy,
        "Ez_xy_nf": Enfz,
        "zmonitor1": zmonitor1,
        "x_xy_nf1": x1,
        "y_xy_nf1": y1,
        "Ex_xy_nf1": Enfx1,
        "Ey_xy_nf1": Enfy1,
        "Ez_xy_nf1": Enfz1,
        "x_ff": x_val,
        "y_ff": y_val,
        "z_ff": z_val,
        "E2_x_ff": E2_x,
        "E2_y_ff": E2_y,
        "E2_z_ff": E2_z,
        "E2_xz_ff": E2_xz,
        "E2_xy_ff": E2_xy
    }, do_compression=True)
    print(f"LUT data saved to {lens_result}")