import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from metalens_portfolio.core import Metalens,grid_and_aperture,hyperbolic_phase,select_elements,make_pupil,scalar_focal_plane

def test_phase_is_radial_and_center_zero():
    l=Metalens();x,y,a=grid_and_aperture(l);p=hyperbolic_phase(x,y,l)
    mid=len(x)//2
    assert np.isclose(p[mid,mid],0)
    assert np.isclose(p[mid,mid+2],p[mid+2,mid])

def test_phase_wrap_nearest():
    r,ph,t,e=select_elements(np.array([2*np.pi-.01]),np.array([.02,np.pi]),np.array([80e-9,150e-9]))
    assert np.isclose(r[0],80e-9)
    assert abs(e[0])<.04

def test_psf_peak_at_center():
    l=Metalens();xx,yy,ap=grid_and_aperture(l)
    fx,fy,psf=scalar_focal_plane(ap.astype(complex),l.pitch_m,l.wavelength_m,l.focal_length_m)
    assert np.unravel_index(psf.argmax(),psf.shape)==tuple(np.array(psf.shape)//2)

def test_pupil():
    mask=np.array([[True,False],[True,True]])
    p=make_pupil(mask,np.zeros((2,2)),np.zeros((2,2)),np.ones((2,2)))
    assert p[0,1]==0 and p[1,1]==1


def test_ideal_airy_fwhm():
    from metalens_portfolio.core import approximate_fwhm
    l = Metalens()
    _, _, ap = grid_and_aperture(l)
    fx, _, psf = scalar_focal_plane(ap, l.pitch_m, l.wavelength_m, l.focal_length_m, oversample=8)
    width = approximate_fwhm(fx, psf[psf.shape[0]//2])
    expected = 1.028*l.wavelength_m*l.focal_length_m/(2*l.radius_m)
    assert abs(width/expected-1) < 0.04
