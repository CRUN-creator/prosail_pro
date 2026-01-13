#!/usr/bin/env python
"""PROSAIL-PRO: PROSPECT-PRO + SAIL radiative transfer models

This package combines:
- PROSPECT-PRO: Leaf optical properties model with protein and CBC parameters
- SAIL: Canopy radiative transfer model (FourSAIL)

PROSPECT-PRO decomposes leaf dry matter (LMA) into:
- Protein content (Prot) in g/cm²
- Carbon-based constituents (CBC) in g/cm²

Relationship: LMA = Prot + CBC

Reference:
Feret, J.-B., et al. (2020). PROSPECT-PRO for estimating content of 
nitrogen-containing leaf proteins and other carbon-based constituents.
Remote Sensing of Environment, 252, 112176.
https://doi.org/10/gh5jqg
"""

__author__ = "J-B Feret (PROSPECT-PRO), J Gomez-Dans (PROSAIL), Python conversion by CRUN"
__version__ = "1.0.0"
__license__ = "GPLv3"

from .spectral_library import get_spectra
spectral_lib = get_spectra()

from .prospect_pro import run_prospect_pro
from .sail_model import run_sail, run_prosail_pro, run_thermal_sail
