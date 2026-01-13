# **PROSAIL-PRO: Python Interface**

**PROSAIL-PRO** is a Python package that combines the **PROSPECT-PRO** leaf optical properties model with the **SAIL** canopy radiative transfer model.

## **Overview**

**PROSPECT-PRO** is the latest evolution of the PROSPECT model family. It extends PROSPECT-D by decomposing leaf dry matter (LMA) into two biologically meaningful components:

* **Proteins (Prot)**  
* **Carbon-based constituents (CBC)** (cellulose, lignin, hemicellulose, sugars, starch, etc.)

The relationship is defined as:

$$
LMA (Cm) = Prot + CBC
$$

This decomposition enables accurate estimation of Nitrogen-related content from leaf optical properties.

## **Features**

* **Dual-Mode Operation**:  
  * **PROSPECT-PRO Mode**: Uses Prot and CBC to simulate leaf properties (Recommended).  
  * **PROSPECT-D Mode**: Uses Cm (LMA) for backward compatibility with legacy research.  
* **Smart Logic**: Automatically detects which mode to use based on input parameters and provides helpful warnings if conflicts arise.  
* **Full PROSAIL Integration**: Seamlessly couples with the 4SAIL canopy model to simulate top-of-canopy reflectance.

## **Installation**

Prerequisites: numpy, scipy

## Install from source  

```
cd prosail_pro
pip install .
```
## Install directly from GitHub

```
pip install git+https://github.com/CRUN-creator/prosail_pro.git
```

## **Usage Guide**

### **1. Running PROSAIL-PRO (Canopy Level)**

To simulate canopy reflectance, use the run_prosail_pro function.

#### **Mode A: PROSPECT-PRO (Standard)**

**Set cm=0** and provide values for prot and cbc.

```
import prosail_pro  
import matplotlib.pyplot as plt

# Simulate canopy reflectance using Proteins and CBC  
rho = prosail_pro.run_prosail_pro(  
    n=1.5,   
    cab=40.0,   
    car=8.0,   
    ant=0.0,   
    cbrown=0.0,     # Note: Parameter name is 'cbrown', not 'brown'  
    cw=0.01,   
    cm=0.0,         # Set Cm to 0 to use Prot + CBC  
    prot=0.001,     # Protein content (g/cm2)  
    cbc=0.009,      # Carbon-based constituents (g/cm2)  
    lai=3.0,   
    lidfa=-0.35,   
    hspot=0.01,   
    tts=30.0,   
    tto=10.0,   
    psi=0.0,   
    rsoil=1.0,   
    psoil=1.0  
)

plt.plot(range(400, 2501), rho)  
plt.title("PROSAIL-PRO Simulation")  
plt.show()
```

#### **Mode B: PROSPECT-D (Legacy/Compatibility)**

**Set prot=0 and cbc=0**, and provide a value for cm.

```
import prosail_pro  
import matplotlib.pyplot as plt

# Simulate using traditional LMA (Cm)  
rho_legacy = prosail_pro.run_prosail_pro(  
    n=1.5, cab=40.0, car=8.0, ant=0.0, cbrown=0.0, cw=0.01,  
    cm=0.01,        # Dry matter content (g/cm2)  
    prot=0.0,       # Set to 0  
    cbc=0.0,        # Set to 0  
    lai=3.0, lidfa=-0.35, hspot=0.01,   
    tts=30.0, tto=10.0, psi=0.0,   
    rsoil=1.0, psoil=1.0  
)

plt.plot(range(400, 2501), rho_legacy)  
plt.title("PROSAIL-D Simulation")  
plt.show()
```

### **2. Running PROSPECT-PRO (Leaf Level Only)**

If you only need leaf reflectance and transmittance:

```
from prosail_pro import run_prospect_pro

# Returns: wavelengths (400-2500nm), Reflectance, Transmittance  
wv, refl, trans = run_prospect_pro(  
    n=1.5, cab=40, car=8, ant=0, cbrown=0,   
    cw=0.01, cm=0, prot=0.001, cbc=0.009  
)
```

### The soil model

The soil model is a fairly simple linear mixture model, where two spectra are mixed and then a brightness term added:

```
rho_soil = rsoil*(psoil*soil_spectrum1+(1-psoil)*soil_spectrum2)
```

The idea is that one of the spectra is a dry soil and the other a wet soil, so soil moisture is then contorlled by `psoil`. `rsoil` is just a brightness scaling term.

## Parameter Reference

### **Leaf Parameters (PROSPECT-PRO)**

| Parameter | Unit | Description | Typical Range |
| :---- | :---- | :---- | ----- |
| **n** | - | Leaf structure parameter | 0.8 - 2.5 |
| **cab** | µg/cm² | Chlorophyll a+b content | 0 - 80 |
| **car** | µg/cm² | Carotenoids content | 0 - 20 |
| **ant** | nmol/cm² | Anthocyanin content | 0 - 40 |
| **cbrown** | - | Brown pigments content (Arbitrary units) | 0 - 1 |
| **cw** | g/cm² | Equivalent water thickness | 0 - 200 |
| **cm** | g/cm² | Dry matter content (LMA). **Set to 0 for PRO mode.** | 0 - 200 |
| **prot** | g/cm² | Protein content. | **0 - 200** |
| **cbc** | g/cm² | Carbon-based constituents. | **0 - 200** |

### **Canopy Parameters (SAIL)**

| Parameter | Unit | Description | Typical Range |
| :---- | :---- | :---- | ----- |
| **lai** | - | Leaf Area Index | 0 - 10 |
| **lidfa** | - | Leaf Angle Distribution | - |
| **lidfb** | - | Leaf Angle Distribution             | - |
| **hspot** | - | Hotspot parameter                   | - |
| **tts** | deg | Solar zenith angle | 0 - 90 |
| **tto** | deg | Observer zenith angle | 0 - 90 |
| **psi** | deg | Relative azimuth angle | 0 - 360 |
| **rsoil** | - | Soil brightness factor (0-1) | 0 - 1 |
| **psoil** | - | Soil moisture factor (0=Wet, 1=Dry) | 0,1 |

### Specifying the leaf angle distribution

The parameter `typelidf` regulates the leaf angle distribution family being used. The following options are understood：

- `typelidf = 1`: use the two parameter LAD parameterisation, where `a` and `b` control the average leaf slope and the distribution bimodality, respectively. Typical distributions are given by the following parameter choices:

| LIDF type    | LIDFa | LIDFb |
| :----------- | :---- | :---- |
| Planophile   | 1     | 0     |
| Erectophile  | -1    | 0     |
| Plagiophile  | 0     | -1    |
| Extremophile | 0     | 1     |
| Spherical    | -0.35 | -0.15 |
| Uniform      | 0     | 0     |

- `typelidf = 2` Ellipsoidal distribution, where `LIDFa` parameter stands for mean leaf angle (0 degrees is planophile, 90 degrees is erectophile). `LIDFb` parameter is ignored.

## **Logic & Warnings**

The model includes built-in logic to handle the relationship between Cm, Prot, and CBC:

1. **Conflict Detected**: If you provide Cm > 0 AND (Prot > 0 OR CBC > 0), the model will **raise a warning** and automatically set Cm = 0 to prioritize the detailed protein/CBC decomposition.  
2. **Info**: If you correctly set Cm=0 (PRO mode) or Prot=0, CBC=0 (D mode), the model will run silently or provide a concise info message depending on configuration.

## **References**

- **Original PROSPECT-PRO Matlab Package**: https://gitlab.com/jbferet/prospect_pro_matlab

- **Original PROSAIL Python Package**: https://github.com/jgomezdans/prosail

## Uninstall

Uninstall the editable version:

```
pip uninstall prosail_pro
```

## Authors

- **CRUN** (Package Maintainer)
- **PROSPECT-PRO (MATLAB)**: [Jean-Baptiste Feret](https://gitlab.com/jbferet)
- **PROSAIL (Python)**: [José Gómez-Dans](https://github.com/jgomezdans)

## Acknowledgments

- Original PROSPECT-PRO MATLAB code by Jean-Baptiste Feret
- Original PROSAIL Python package by José Gómez-Dans
