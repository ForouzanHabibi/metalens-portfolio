"""Run a license-free demonstration using SYNTHETIC, non-FDTD phase-library data."""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from metalens_portfolio.core import (Metalens, grid_and_aperture, hyperbolic_phase,
                                   select_elements, make_pupil, scalar_focal_plane,
                                   approximate_fwhm)

OUT = Path(__file__).resolve().parents[1] / 'figures'
OUT.mkdir(exist_ok=True)
lens = Metalens()
xx, yy, aperture = grid_and_aperture(lens)
target = hyperbolic_phase(xx, yy, lens)
# Synthetic library for demonstrating the software only; NOT physical FDTD predictions.
radii = np.linspace(65e-9, 220e-9, 24)
phase = np.linspace(0, 2*np.pi*0.97, len(radii))
transmission = 0.72 + 0.18 * np.sin(np.linspace(0, np.pi, len(radii)))**2
chosen_radii, chosen_phase, chosen_T, err = select_elements(target, phase, radii, transmission)
ideal = aperture.astype(complex)
realized = make_pupil(aperture, target, chosen_phase, chosen_T)
fx, fy, psf_ideal = scalar_focal_plane(ideal, lens.pitch_m, lens.wavelength_m, lens.focal_length_m)
_, _, psf = scalar_focal_plane(realized, lens.pitch_m, lens.wavelength_m, lens.focal_length_m)
cy, cx = np.array(psf.shape)//2
width = approximate_fwhm(fx, psf[cy, :])
print(f'Design NA: {lens.na:.3f}; sampling bound lambda/(2p): {lens.phase_sampling_na_limit:.3f}')
print(f'Synthetic-library phase RMS error: {np.sqrt(np.mean(err[aperture]**2)):.3f} rad')
print(f'Scalar FWHM (synthetic design): {width*1e6:.3f} micrometers')
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(radii*1e9, phase, 'o-', label='Synthetic phase')
ax.set(xlabel='Pillar radius (nm)', ylabel='Phase (rad)', title='Illustrative phase library — not FDTD')
ax.legend();fig.tight_layout();fig.savefig(OUT/'synthetic_phase_library.png',dpi=170);plt.close(fig)
fig,axs = plt.subplots(1,2,figsize=(9,4))
maskrad=np.where(aperture, chosen_radii*1e9, np.nan)
im=axs[0].imshow(maskrad, origin='lower', extent=[xx.min()*1e6,xx.max()*1e6,yy.min()*1e6,yy.max()*1e6])
fig.colorbar(im,ax=axs[0],label='Radius (nm)');axs[0].set_title('Designed pillar radii')
axs[1].imshow(np.where(aperture,err,np.nan),origin='lower',vmin=-np.pi,vmax=np.pi)
axs[1].set_title('Wrapped phase mismatch (rad)')
fig.tight_layout();fig.savefig(OUT/'synthetic_metalens_layout.png',dpi=170);plt.close(fig)
fig,ax=plt.subplots(figsize=(6,4));extent=[fx[0]*1e6,fx[-1]*1e6,fy[0]*1e6,fy[-1]*1e6]
ax.imshow(psf,origin='lower',extent=extent,vmin=0,vmax=0.35)
ax.set(xlabel='x (um)',ylabel='y (um)',title='Illustrative normalized scalar PSF')
ax.set_xlim(-7,7);ax.set_ylim(-7,7)
fig.tight_layout();fig.savefig(OUT/'synthetic_psf.png',dpi=170);plt.close(fig)
np.savetxt(OUT/'synthetic_phase_library.csv',np.column_stack([radii,phase,transmission]),delimiter=',',header='radius_m,phase_rad,power_transmission',comments='')
