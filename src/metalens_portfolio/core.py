"""Independent numerical utilities for phase-library metalens prototyping.

The functions in this module are new example implementations and do not
contain or require unpublished collaborator scripts or proprietary data.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Metalens:
    wavelength_m: float = 1.55e-6
    focal_length_m: float = 40e-6
    radius_m: float = 12e-6
    pitch_m: float = 0.5e-6

    def __post_init__(self):
        if min(self.wavelength_m, self.focal_length_m, self.radius_m, self.pitch_m) <= 0:
            raise ValueError('All geometry and wavelength values must be positive')

    @property
    def na(self):
        return self.radius_m / np.hypot(self.radius_m, self.focal_length_m)

    @property
    def phase_sampling_na_limit(self):
        return self.wavelength_m / (2 * self.pitch_m)


def grid_and_aperture(lens: Metalens):
    count = int(np.ceil(lens.radius_m / lens.pitch_m))
    coords = np.arange(-count, count + 1) * lens.pitch_m
    xx, yy = np.meshgrid(coords, coords, indexing='xy')
    aperture = xx**2 + yy**2 <= lens.radius_m**2
    return xx, yy, aperture


def hyperbolic_phase(xx, yy, lens: Metalens):
    """Ideal on-axis focusing phase, modulo 2 pi."""
    path_delta = np.sqrt(xx**2 + yy**2 + lens.focal_length_m**2) - lens.focal_length_m
    return np.mod(-2 * np.pi * path_delta / lens.wavelength_m, 2 * np.pi)


def select_elements(target_phase, library_phase, library_radius, library_transmission=None):
    """Select meta-atoms by wrapped phase error (no optimization claimed).

    If supplied, library_transmission is carried through as *power* transmission.
    """
    phi = np.asarray(library_phase, dtype=float).reshape(-1)
    rad = np.asarray(library_radius, dtype=float).reshape(-1)
    if len(phi) < 2 or len(phi) != len(rad):
        raise ValueError('Provide matching phase and radius arrays with >=2 entries')
    if not np.all(np.isfinite(phi)) or not np.all(np.isfinite(rad)) or np.any(rad <= 0):
        raise ValueError('Library must have finite phase and positive radii')
    if library_transmission is None:
        transmission = np.ones_like(phi)
    else:
        transmission = np.asarray(library_transmission, dtype=float).reshape(-1)
        if transmission.shape != phi.shape or np.any((transmission < 0) | (transmission > 1)):
            raise ValueError('Power transmission must be within [0, 1] and match phases')
    phase_error = np.angle(np.exp(1j * (np.asarray(target_phase)[..., None] - phi)))
    ids = np.argmin(np.abs(phase_error), axis=-1)
    return rad[ids], phi[ids], transmission[ids], np.take_along_axis(phase_error, ids[..., None], axis=-1)[..., 0]


def scalar_focal_plane(pupil, pitch_m, wavelength_m, focal_length_m, oversample=4):
    """Fraunhofer scalar intensity for a pupil sampled on a Cartesian grid.

    An illustrative thin-mask model, *not* an FDTD simulation or coronagraph.
    Returned intensity is normalized to peak 1 for easy visualization.
    """
    a = np.asarray(pupil, dtype=complex)
    if a.ndim != 2 or min(a.shape) < 3 or oversample < 1:
        raise ValueError('Expected a 2D pupil and positive oversampling')
    shape = tuple(int(d * oversample) for d in a.shape)
    # Pad the contiguous pupil before transforming. For intensity only,
    # translating the pupil changes Fourier phase but not intensity.
    field = np.fft.fftshift(np.fft.fft2(a, s=shape))
    intensity = np.abs(field)**2
    intensity /= max(float(intensity.max()), np.finfo(float).tiny)
    fx = np.fft.fftshift(np.fft.fftfreq(shape[1], d=pitch_m))
    fy = np.fft.fftshift(np.fft.fftfreq(shape[0], d=pitch_m))
    return fx * wavelength_m * focal_length_m, fy * wavelength_m * focal_length_m, intensity


def make_pupil(aperture, target_phase, realized_phase, power_transmission, phase_error=None):
    """Return residual-phase pupil for focus-plane comparison.

    Ideal focusing phase is analytically removed: it represents propagation
    to the ideal focal plane in the scalar approximation.
    """
    residual = np.angle(np.exp(1j * (realized_phase - target_phase)))
    if phase_error is not None:
        residual = np.asarray(phase_error)
    return np.asarray(aperture, dtype=float) * np.sqrt(power_transmission) * np.exp(1j * residual)


def approximate_fwhm(coords, intensity_line):
    """Approximate FWHM of central lobe using interpolated half-max crossings."""
    x = np.asarray(coords, dtype=float)
    y = np.asarray(intensity_line, dtype=float)
    center = int(np.argmax(y))
    half = y[center] / 2
    left, right = center, center
    while left > 0 and y[left] >= half:
        left -= 1
    while right < len(y)-1 and y[right] >= half:
        right += 1
    if left == 0 or right == len(y)-1:
        return np.nan
    xl = np.interp(half, [y[left], y[left+1]], [x[left], x[left+1]])
    xr = np.interp(half, [y[right], y[right-1]], [x[right], x[right-1]])
    return float(xr-xl)
