#!/usr/bin/env python
"""Spectral libraries for PROSAIL-PRO

This module loads spectral data required by SAIL model:
- Soil reflectance spectra (dry and wet)
- Light spectra (direct and diffuse)

Note: PROSPECT-PRO spectral data is loaded separately in prospect_pro.py
"""
import pkgutil
import os
from collections import namedtuple
from io import BytesIO
import numpy as np

Spectra = namedtuple('Spectra', 'prospect_pro soil light')
ProspectProSpectra = namedtuple('ProspectProSpectra', 
                                'nr kab kcar kant kbrown kw km kprot kcbc')
SoilSpectra = namedtuple("SoilSpectra", "rsoil1 rsoil2")
LightSpectra = namedtuple("LightSpectra", "es ed")

def _load_data_safe(filename):
    """Helper function: Attempt to load from the package, fall back to loading from the local directory if failed"""
    try:
        data = pkgutil.get_data('prosail_pro', filename)
    except ImportError:
        data = None
    
    if data is None:
        # 本地回退机制
        current_dir = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(current_dir, filename)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Cannot find {filename} in package or {current_dir}")
        return np.loadtxt(file_path, unpack=True)
    
    return np.loadtxt(BytesIO(data), unpack=True)

def get_spectra():
    """Reads the spectral information for SAIL model"""
    # PROSPECT-PRO
    p_data = _load_data_safe('prospect_pro_spectra.txt')
    
    _, nr, kab, kcar, kant, kbrown, kw, km, kprot, kcbc = p_data
    prospect_pro_spectra = ProspectProSpectra(nr, kab, kcar, kant, kbrown, kw, km, kprot, kcbc)

    # SOIL
    rsoil1, rsoil2 = _load_data_safe('soil_reflectance.txt')
    soil_spectra = SoilSpectra(rsoil1, rsoil2)    
    
    # LIGHT
    es, ed = _load_data_safe('light_spectra.txt')
    light_spectra = LightSpectra(es, ed)
    
    return Spectra(prospect_pro_spectra, soil_spectra, light_spectra)
