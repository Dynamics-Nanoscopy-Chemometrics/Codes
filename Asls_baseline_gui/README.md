# ASLS Baseline GUI: Asymmetric Least-Squares Baseline Correction

Interactive MATLAB GUI for **global and local** baseline correction using Asymmetric Least-Squares (ASLS) smoothing.

**Reference**:
Eilers, P. H. C., & Boelens, H. F. M. (2005). Baseline correction with asymmetric least squares smoothing. *Leiden University Medical Centre Report*, 1(1), 5.

---

## Overview

The **ASLS Baseline GUI** provides an interactive tool for removing baseline drift from spectroscopic data using the Asymmetric Least-Squares (ASLS) method. This implementation extends the classic ASLS approach by supporting both **global** and **local** baseline corrections, allowing users to apply different correction parameters to specific spectral regions.

The baseline is computed by:
1. Estimating a **global baseline** using ASLS across the entire spectrum
2. Computing **local corrections** on user-specified intervals with custom parameters
3. Merging and blending these corrections to form a single smooth baseline
4. Optionally applying **final smoothing** to the merged baseline

This approach is particularly useful for complex spectral data where baseline characteristics vary across different wavelength regions.


---

## Algorithm

### Asymmetric Least-Squares (ASLS)

The ASLS method estimates a baseline **z** from signal **y** by solving:

**min ||W(y - z)||² + λ||Dz||²**

where:
- **W**: Diagonal weight matrix (updated iteratively)
- **D**: Second-difference matrix (penalizes curvature)
- **λ**: Smoothness parameter (larger = smoother baseline)
- **p**: Asymmetry parameter (controls weighting: 0 < p < 1)

### Iterative Weighting

At each iteration:
1. Solve for baseline **z** using current weights
2. Update weights asymmetrically:
   - **w[i] = p** if y[i] > z[i] (data above baseline)
   - **w[i] = 1-p** if y[i] < z[i] (data below baseline)

This asymmetric weighting allows the baseline to pass below peaks while fitting the actual baseline.

---

## Installation

### Prerequisites
- MATLAB R2016a or later
- No additional toolboxes required

### Setup

1. Add the folder to your MATLAB path:

```matlab
addpath('path/to/Codes/Asls_baseline_gui');
```

2. Verify installation:

```matlab
which AsLSLocalGlobalGUI
which AslsLocalParts
```

### Dependencies

All required functions are included in this folder:
- `AsLSLocalGlobalGUI.m` — Interactive GUI
- `AslsLocalParts.m` — Standalone batch function
- `AslsBaseSingle.m` — Core ASLS solver
- `LocalBaselineOnly.m` — Local interval correction

---

## Function Reference

### `AsLSLocalGlobalGUI`

```matlab
AsLSLocalGlobalGUI
```

Launches the interactive GUI for ASLS baseline correction with global and local corrections.

---

## Parameter Selection Guidelines

### Lambda (Smoothness)
- **Small λ (1e3 - 1e5)**: Follows signal closely, captures fine baseline variations
- **Medium λ (1e5 - 1e7)**: Balanced smoothness, most common choice
- **Large λ (1e7 - 1e9)**: Very smooth baseline, removes broad features

### p (Asymmetry)
- **Small p (0.001 - 0.01)**: Baseline passes well below peaks (recommended)
- **Medium p (0.01 - 0.1)**: More symmetric weighting
- **Large p (0.1 - 0.5)**: Approaches symmetric smoothing

### General Recommendations
1. Start with **λ = 1e6, p = 0.001** for global parameters
2. Use **Preview Baseline** to assess fit on mean spectrum
3. Add local intervals for regions with distinct baseline features
4. Use **smaller λ** in local intervals for more responsive correction
5. Enable **final smoothing** (λ = 1e5) for seamless transitions

---

## License

Released under the **MIT License**.

---

## Authors

- **Yesid Roman Gómez**
- **Date Created**: October 16, 2025
- **Reviewed by**: Lovelace's Square

---

## Changelog

- **v1.0.0 (2025-10-16)**:
  Initial release with GUI and standalone function for global + local ASLS baseline correction