#!/usr/bin/env python
"""SAIL radiative transfer model integrated with PROSPECT-PRO

This module combines PROSPECT-PRO leaf optical properties model with
the SAIL canopy radiative transfer model.
"""
import numpy as np

from .prospect_pro import run_prospect_pro
from .FourSAIL import foursail
from .spectral_library import get_spectra

_SPECTRA_CACHE = get_spectra()

def run_prosail_pro(n, cab, car, ant, cbrown, cw, cm, prot, cbc,
                    lai, lidfa, hspot, tts, tto, psi,
                    alpha=40., typelidf=2, lidfb=0., factor="SDR",
                    rsoil0=None, rsoil=None, psoil=None,
                    soil_spectrum1=None, soil_spectrum2=None):
    """Run PROSPECT-PRO + SAIL radiative transfer models
    
    This function combines PROSPECT-PRO for leaf optical properties with
    SAIL for canopy radiative transfer.
    
    Parameters
    ----------
    n : float
        Leaf structure parameter
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
    lai : float
        Leaf area index
    lidfa : float
        Leaf angle distribution parameter a. If typelidf=2, average leaf inclination angle
    hspot : float
        Hotspot parameter
    tts : float
        Solar zenith angle (degrees)
    tto : float
        Sensor zenith angle (degrees)
    psi : float
        Relative sensor-solar azimuth angle (degrees, saa - vaa)
    alpha : float, optional
        Angle for surface scattering (default: 40 degrees)
    typelidf : int, optional
        Type of leaf angle distribution function (default: 2)
    lidfb : float, optional
        Leaf angle distribution parameter b (ignored if typelidf=2)
    factor : str, optional
        Reflectance factor to return (default: "SDR")
        - "SDR": directional reflectance factor
        - "BHR": bi-hemispherical reflectance factor
        - "DHR": directional-hemispherical reflectance factor
        - "HDR": hemispherical-directional reflectance factor
        - "ALL": all of the above
        - "ALLALL": all terms calculated by SAIL
    rsoil0 : array-like, optional
        Soil reflectance spectrum (2101 elements)
    rsoil : float, optional
        Soil brightness scalar
    psoil : float, optional
        Soil moisture scalar
    soil_spectrum1 : array-like, optional
        First soil spectrum component (2101 elements)
    soil_spectrum2 : array-like, optional
        Second soil spectrum component (2101 elements)
        
    Returns
    -------
    array or list of arrays
        Reflectance factor(s) between 400 and 2500 nm, depending on 'factor' parameter
        
    Notes
    -----
    LMA (Leaf Mass per Area) = Prot + CBC
    
    The soil model is a linear mixture:
        rho_soil = rsoil * (psoil * soil_spectrum1 + (1-psoil) * soil_spectrum2)
    
    By default, soil_spectrum1 is dry soil and soil_spectrum2 is wet soil.
    
    Examples
    --------
    >>> import prosail_pro
    >>> canopy_refl = prosail_pro.run_prosail_pro(
    ...     n=1.5, cab=40.0, car=8.0, ant=0.5, cbrown=0.0,
    ...     cw=0.01, cm=0.1, prot=0.001, cbc=0.009,
    ...     lai=3.0, lidfa=57.0, hspot=0.01,
    ...     tts=30.0, tto=10.0, psi=0.0,
    ...     rsoil=0.15, psoil=1.0
    ... )
    """
    
    factor = factor.upper()
    if factor not in ["SDR", "BHR", "DHR", "HDR", "ALL", "ALLALL"]:
        raise ValueError(
            "'factor' must be one of SDR, BHR, DHR, HDR, ALL or ALLALL"
        )
    
    # Set default soil spectra
    if soil_spectrum1 is not None:
        assert len(soil_spectrum1) == 2101, "soil_spectrum1 must have 2101 elements"
    else:
        soil_spectrum1 = _SPECTRA_CACHE.soil.rsoil1
    
    if soil_spectrum2 is not None:
        assert len(soil_spectrum2) == 2101, "soil_spectrum2 must have 2101 elements"
    else:
        soil_spectrum2 = _SPECTRA_CACHE.soil.rsoil2
    
    # Calculate soil reflectance
    if rsoil0 is None:
        if (rsoil is None) or (psoil is None):
            raise ValueError(
                "If rsoil0 isn't defined, then rsoil and psoil need to be defined!"
            )
        rsoil0 = rsoil * (psoil * soil_spectrum1 + (1. - psoil) * soil_spectrum2)
    
    # Run PROSPECT-PRO to get leaf optical properties
    wv, refl, trans = run_prospect_pro(
        n=n, cab=cab, car=car, ant=ant, cbrown=cbrown, 
        cw=cw, cm=cm, prot=prot, cbc=cbc, alpha=alpha
    )
    
    # Run SAIL model
    [tss, too, tsstoo, rdd, tdd, rsd, tsd, rdo, tdo,
     rso, rsos, rsod, rddt, rsdt, rdot, rsodt, rsost, rsot,
     gammasdf, gammasdb, gammaso] = foursail(refl, trans,
                                              lidfa, lidfb, typelidf,
                                              lai, hspot,
                                              tts, tto, psi, rsoil0)
    
    # Return requested reflectance factor
    if factor == "SDR":
        return rsot
    elif factor == "BHR":
        return rddt
    elif factor == "DHR":
        return rsdt
    elif factor == "HDR":
        return rdot
    elif factor == "ALL":
        return [rsot, rddt, rsdt, rdot]
    elif factor == "ALLALL":
        return [tss, too, tsstoo, rdd, tdd, rsd, tsd, rdo, tdo,
                rso, rsos, rsod, rddt, rsdt, rdot, rsodt, rsost, rsot,
                gammasdf, gammasdb, gammaso]


def run_sail(refl, trans, lai, lidfa, hspot, tts, tto, psi,
             typelidf=2, lidfb=0., factor="SDR",
             rsoil0=None, rsoil=None, psoil=None,
             soil_spectrum1=None, soil_spectrum2=None):
    """Run SAIL radiative transfer model with custom leaf spectra
    
    This function allows you to run SAIL with pre-calculated leaf reflectance
    and transmittance from PROSPECT-PRO or any other source.
    
    Parameters
    ----------
    refl : array-like
        Leaf reflectance (2101 elements, 400-2500 nm)
    trans : array-like
        Leaf transmittance (2101 elements, 400-2500 nm)
    lai : float
        Leaf area index
    lidfa : float
        Leaf angle distribution parameter a
    hspot : float
        Hotspot parameter
    tts : float
        Solar zenith angle (degrees)
    tto : float
        Sensor zenith angle (degrees)
    psi : float
        Relative sensor-solar azimuth angle (degrees)
    typelidf : int, optional
        Type of leaf angle distribution function (default: 2)
    lidfb : float, optional
        Leaf angle distribution parameter b
    factor : str, optional
        Reflectance factor to return (default: "SDR")
    rsoil0 : array-like, optional
        Soil reflectance spectrum
    rsoil : float, optional
        Soil brightness scalar
    psoil : float, optional
        Soil moisture scalar
    soil_spectrum1 : array-like, optional
        First soil spectrum component
    soil_spectrum2 : array-like, optional
        Second soil spectrum component
        
    Returns
    -------
    array or list of arrays
        Reflectance factor(s) between 400 and 2500 nm
        
    Examples
    --------
    >>> import prosail_pro
    >>> # First get leaf optical properties
    >>> wv, refl, trans = prosail_pro.run_prospect_pro(
    ...     n=1.5, cab=40.0, car=8.0, ant=0.5, brown=0.0,
    ...     cw=0.01, prot=0.001, cbc=0.009
    ... )
    >>> # Then run SAIL
    >>> canopy_refl = prosail_pro.run_sail(
    ...     refl, trans, lai=3.0, lidfa=57.0, hspot=0.01,
    ...     tts=30.0, tto=10.0, psi=0.0,
    ...     rsoil=0.15, psoil=1.0
    ... )
    """
    from . import spectral_lib
    
    factor = factor.upper()
    if factor not in ["SDR", "BHR", "DHR", "HDR", "ALL", "ALLALL"]:
        raise ValueError(
            "'factor' must be one of SDR, BHR, DHR, HDR, ALL or ALLALL"
        )
    
    # Set default soil spectra
    if soil_spectrum1 is not None:
        assert (len(soil_spectrum1) == 2101)
    else:
        soil_spectrum1 = spectral_lib.soil.rsoil1
    
    if soil_spectrum2 is not None:
        assert (len(soil_spectrum2) == 2101)
    else:
        soil_spectrum2 = spectral_lib.soil.rsoil2
    
    # Calculate soil reflectance
    if rsoil0 is None:
        if (rsoil is None) or (psoil is None):
            raise ValueError(
                "If rsoil0 isn't defined, then rsoil and psoil need to be defined!"
            )
        rsoil0 = rsoil * (psoil * soil_spectrum1 + (1. - psoil) * soil_spectrum2)
    
    # Run SAIL model
    [tss, too, tsstoo, rdd, tdd, rsd, tsd, rdo, tdo,
     rso, rsos, rsod, rddt, rsdt, rdot, rsodt, rsost, rsot,
     gammasdf, gammasdb, gammaso] = foursail(refl, trans,
                                              lidfa, lidfb, typelidf,
                                              lai, hspot,
                                              tts, tto, psi, rsoil0)
    
    # Return requested reflectance factor
    if factor == "SDR":
        return rsot
    elif factor == "BHR":
        return rddt
    elif factor == "DHR":
        return rsdt
    elif factor == "HDR":
        return rdot
    elif factor == "ALL":
        return [rsot, rddt, rsdt, rdot]
    elif factor == "ALLALL":
        return [tss, too, tsstoo, rdd, tdd, rsd, tsd, rdo, tdo,
                rso, rsos, rsod, rddt, rsdt, rdot, rsodt, rsost, rsot,
                gammasdf, gammasdb, gammaso]

def run_thermal_sail(lam,  
                     tveg, tsoil, tveg_sunlit, tsoil_sunlit, t_atm, 
                     lai, lidfa, hspot,  
                     tts, tto, psi, rsoil=None,
                     refl=None, emv=None, ems=None,
                     typelidf=2, lidfb=0):
    c1 = 3.741856E-16
    c2 = 14388.0
    # Calculate the thermal emission from the different
    # components using Planck's Law
    top = (1.0e-6)*c1*(lam*1e-6)**(-5.)
    Hc = top / ( np.exp ( c2/(lam*tveg))-1.)         # Shade leaves
    Hh = top / ( np.exp ( c2/(lam*tveg_sunlit))-1.)  # Sunlit leaves
    Hd = top / ( np.exp ( c2/(lam*tsoil))-1.)        # shade soil 
    Hs = top / ( np.exp ( c2/(lam*tsoil_sunlit))-1.) # Sunlit soil
    Hsky = top / ( np.exp ( c2/(lam*t_atm))-1.)      # Sky emission
    
    
    # Emissivity calculations
    if refl is not None and emv is None:
        emv = 1. - refl # Assuming absorption is 1
    
    if rsoil is not None and ems is None:
        ems = 1. - rsoil
    
    if rsoil is None and ems is not None:
        rsoil = 1. - ems
    if refl is None and emv is not None:
        refl = 1. - emv
    
    [tss, too, tsstoo, rdd, tdd, rsd, tsd, rdo, tdo,
         rso, rsos, rsod, rddt, rsdt, rdot, rsodt, rsost, rsot,
         gammasdf, gammasdb, gammaso] = foursail (refl, np.zeros_like(refl),  
                                                  lidfa, lidfb, typelidf, 
                                                  lai, hspot, 
                                                  tts, tto, psi, rsoil)
    
    gammad = 1.0 - rdd - tdd
    gammao = 1.0 - rdo - tdo - too

    tso = tss*too+tss*(tdo+rsoil*rdd*too)/(1.0-rsoil*rdd)
    ttot = (too+tdo)/(1.0-rsoil*rdd)
    gammaot = gammao + ttot*rsoil*gammad
    gammasot = gammaso + ttot*rsoil*gammasdf

    aeev = gammaot
    aees = ttot*ems

    Lw = ( rdot*Hsky + 
            (aeev*Hc + 
            gammasot*emv*(Hh-Hc) + 
            aees*Hd + 
            tso*ems*(Hs-Hd)))/np.pi
    
    dnoem1 = top/(Lw*np.pi)
    Tbright = c2/(lam*np.log(dnoem1+1.0))
    dir_em = 1.0 - rdot 
    return Lw, Tbright, dir_em