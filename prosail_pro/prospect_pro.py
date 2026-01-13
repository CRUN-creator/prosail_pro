#!/usr/bin/env python
"""The PROSPECT-PRO leaf optical properties model
PROSPECT-PRO extends PROSPECT-D by decomposing dry matter into:
- Protein content (Prot) in g/cm²
- Carbon-based constituents (CBC) in g/cm²
Relationship: LMA = Prot + CBC
Reference:
Feret et al. (2020) PROSPECT-PRO for estimating content of
nitrogen-containing leaf proteins and other carbon-based constituents
Remote Sensing of Environment, 252, 112176.
https://doi.org/10.1016/j.rse.2020.112176
"""
import numpy as np
from scipy.special import expi
import os
import warnings
from . import spectral_library
from .spectral_library import get_spectra
_SPECTRA_CACHE = get_spectra().prospect_pro

def run_prospect_pro(n, cab, car, ant, cbrown, cw, cm, prot, cbc, 
                     nr=None, kab=None, kcar=None, kbrown=None, kw=None, 
                     km=None, kprot=None, kcbc=None, kant=None, alpha=40.):
    """
    The PROSPECT-PRO model interface.
    
    Acts as a wrapper to handle spectral data loading automatically,
    mimicking the style of the official run_prospect function.
    """
    
    # Get the default spectral library

    libs = spectral_library.get_spectra().prospect_pro
    
    # Call the core calculation function
    # If the user does not provide a specific spectral coefficient (None), the default value from the library will be used.
    wv, refl, trans = prospect_pro(
        n, cab, car, ant, cbrown, cw, cm, prot, cbc,
        _SPECTRA_CACHE.nr if nr is None else nr,
        _SPECTRA_CACHE.kab if kab is None else kab,
        _SPECTRA_CACHE.kcar if kcar is None else kcar,
        _SPECTRA_CACHE.kant if kant is None else kant,
        _SPECTRA_CACHE.kbrown if kbrown is None else kbrown,
        _SPECTRA_CACHE.kw if kw is None else kw,
        _SPECTRA_CACHE.km if km is None else km,
        _SPECTRA_CACHE.kprot if kprot is None else kprot,
        _SPECTRA_CACHE.kcbc if kcbc is None else kcbc,
        alpha=alpha
    )

    return wv, refl, trans

def calctav(alpha, nr):
    """Calculate transmittance at the interface
    Based on:
    - Stern F. (1964), Transmission of isotropic radiation across an
      interface between two dielectrics, Appl. Opt., 3(1):111-113.
    - Allen W.A. (1973), Transmission of isotropic light across a
      dielectric surface in two and three dimensions, J. Opt. Soc. Am.,
      63(6):664-666.
    Parameters
    ----------
    alpha : float
        Angle in degrees
    nr : array-like
        Refractive index
    Returns
    -------
    array-like
        Transmittance values
    """
    #rd  = pi/180 np.deg2rad
    n2 = nr * nr
    npx = n2 + 1
    nm = n2 - 1
    a = (nr + 1) * (nr + 1) / 2.
    k = -(n2 - 1) * (n2 - 1) / 4.
    sa = np.sin(np.deg2rad(alpha))

    if alpha != 90:
        b1 = np.sqrt((sa*sa - npx/2) * (sa*sa - npx/2) + k)
    else:
        b1 = 0.
    b2 = sa*sa - npx/2
    b = b1 - b2
    b3 = b**3
    a3 = a**3
    ts = (k**2/(6*b3) + k/b - b/2) - (k**2./(6*a3) + k/a - a/2)
    tp1 = -2*n2*(b - a)/(npx**2)
    tp2 = -2*n2*npx*np.log(b/a)/(nm**2)
    tp3 = n2*(1/b - 1/a)/2
    tp4 = 16*n2**2*(n2**2 + 1)*np.log((2*npx*b - nm**2)/(2*npx*a - nm**2))/(npx**3*nm**2)
    tp5 = 16*n2**3*(1./(2*npx*b - nm**2) - 1/(2*npx*a - nm**2))/(npx**3)
    tp = tp1 + tp2 + tp3 + tp4 + tp5
    tav = (ts + tp)/(2*sa**2)

    return tav

def refl_trans_one_layer(alpha, nr, tau):
    """Calculate reflectance and transmittance of one layer
    Based on:
    Allen W.A., Gausman H.W., Richardson A.J., Thomas J.R. (1969),
    Interaction of isotropic light with a compact plant leaf, J. Opt.
    Soc. Am., 59(10):1376-1379.
    Parameters
    ----------
    alpha : float
        Angle in degrees
    nr : array-like
        Refractive index
    tau : array-like
        Transmittance of the layer
    Returns
    -------
    tuple
        (r, t, Ra, Ta, denom) - reflectance and transmittance values
    """
    # Reflectivity and transmissivity at the interface
    talf = calctav(alpha, nr)
    ralf = 1.0 - talf
    t12 = calctav(90, nr)
    r12 = 1. - t12
    t21 = t12 / (nr * nr)
    r21 = 1 - t21

    # Top surface side
    denom = 1. - r21*r21*tau*tau
    Ta = talf*tau*t21/denom
    Ra = ralf + r21*tau*Ta

    # Bottom surface side
    t = t12*tau*t21/denom
    r = r12 + r21*tau*t

    return r, t, Ra, Ta, denom

def prospect_pro(N, cab, car, ant, cbrown, cw, cm, prot, cbc,
                 nr, kab, kcar, kant, kbrown,
                 kw, km, kprot, kcbc, alpha=40.):
    """PROSPECT-PRO model for leaf optical properties
    Parameters
    ----------
    N : float
        Leaf structure parameter (number of layers)
    cab : float
        Chlorophyll a+b content (µg/cm²)
    car : float
        Carotenoids content (µg/cm²)
    ant : float
        Anthocyanin content (nmol/cm² or µg/cm²)
    cbrown : float
        Brown pigments content (arbitrary units)
    cw : float
        Equivalent water thickness (g/cm² or cm)
    cm : float
        Leaf dry matter content (g/cm²)
    prot : float
        Protein content (g/cm²)
    cbc : float
        Carbon-based constituents content (g/cm²)
    nr : array-like, optional
        Refractive index (2101 elements, 400-2500 nm)
    kab : array-like, optional
        Chlorophyll absorption coefficient
    kcar : array-like, optional
        Carotenoid absorption coefficient
    kant : array-like, optional
        Anthocyanin absorption coefficient
    kbrown : array-like, optional
        Brown pigment absorption coefficient
    kw : array-like, optional
        Water absorption coefficient
    km : array-like, optional
        Dry matter absorption coefficient
    kprot : array-like, optional
        Protein absorption coefficient
    kcbc : array-like, optional
        CBC absorption coefficient
    alpha : float, optional
        Angle for surface scattering (default: 40 degrees)  
    Returns
    -------
    tuple
        (wavelengths, reflectance, transmittance)
    Notes
    -----
    LMA (Leaf Mass per Area) = Prot + CBC
    """
    lambdas = np.arange(400, 2501)
    n_lambdas = len(lambdas)
    n_elems_list = [len(spectrum)  for spectrum in
                [nr, kab, kcar, kant, kbrown, kw, km, kprot, kcbc]]
    if not all(n_elems == n_lambdas for n_elems in n_elems_list):
        raise ValueError("Leaf spectra don't have the right shape!")

    if cm > 0 and (cbc > 0 or prot > 0):
        warnings.warn('PROSPECT-PRO expects Leaf Mass Per Area (LMA, Cm parameter) '
                      'to be set to 0, as LMA is decomposed into CBC and Proteins. '
                      'Cm is automatically set to 0 to respect the relationship '
                      'Cm = CBC + Prot.')
        warnings.warn('Please set Prot and CBC to 0 if you want to use PROSPECT-D')
        cm = 0
    elif cm > 0 and cbc == 0 and prot == 0:
        warnings.warn('[Info] You are using PROSPECT-D mode (Calculation based on Cm/LMA).')
    elif cm == 0 and (cbc > 0 or prot > 0):
        warnings.warn('[Info] You are using PROSPECT-PRO mode (Calculation based on Prot + CBC).')


    # Calculate total absorption
    kall = (cab*kab + car*kcar + ant*kant + cbrown*kbrown +
            cw*kw + cm*km + prot*kprot + cbc*kcbc) / N
    
    # Non-conservative scattering (normal case)
    j = kall > 0
    t1 = (1 - kall) * np.exp(-kall)
    t2 = kall**2 * (-expi(-kall))
    tau = np.ones_like(t1)
    tau[j] = t1[j] + t2[j]

    # Reflectance and transmittance of one layer
    r, t, Ra, Ta, denom = refl_trans_one_layer(alpha, nr, tau)

    # Reflectance and transmittance of N layers
    # Stokes equations (Stokes G.G., 1862)
    D = np.sqrt((1 + r + t) * (1 + r - t) * (1. - r + t) * (1. - r - t))
    rq = r * r
    tq = t * t
    a = (1 + rq - tq + D) / (2 * r)
    b = (1 - rq + tq + D) / (2 * t)

    bNm1 = np.power(b, N - 1)
    bN2 = bNm1 * bNm1
    a2 = a * a
    denom = a2 * bN2 - 1
    Rsub = a * (bN2 - 1) / denom
    Tsub = bNm1 * (a2 - 1) / denom

    # Case of zero absorption
    j = r + t >= 1.
    Tsub[j] = t[j] / (t[j] + (1 - t[j]) * (N - 1))
    Rsub[j] = 1 - Tsub[j]

    # Reflectance and transmittance of the leaf: combine top layer with next N-1 layers
    denom = 1 - Rsub * r
    tran = Ta * Tsub / denom
    refl = Ra + Ta * Rsub * t / denom

    return lambdas, refl, tran
