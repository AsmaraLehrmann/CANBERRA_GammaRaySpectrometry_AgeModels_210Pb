# Mass Accumulation Rate (MAR) Calculator

Standalone Python script for computing mass accumulation rates from a 210Pb
age model and a wet/dry weights table. Companion to
`210Pb_AgeModelPlot.ipynb` and `sedimentation_rate.py`.

## What it does

1. Reads an **AgeModel CSV** (the file produced by `210Pb_AgeModelPlot.ipynb`)
   and a **weights CSV** that has wet and dry sediment masses per interval.
2. Computes **dry bulk density (DBD)** for each interval using the indirect
   Sanchez-Cabeza & Ruiz-Fernandez (2012) method:

   $$\rho_{db} = \frac{M_d}{\dfrac{M_w - M_d}{\rho_w} + \dfrac{M_d}{\rho_s}}$$

   with $\rho_w = 1.02$ g/cm³ (seawater) and $\rho_s = 2.65$ g/cm³ (quartz).
3. Computes **effective interval thickness** using the midpoint-to-midpoint
   convention (option b from our discussion): each sample's $\Delta z$ runs
   from the midpoint with its shallower neighbor to the midpoint with its
   deeper neighbor (endpoints clamped to the top of the first interval and
   the base of the last interval). This fills in gaps between non-contiguous
   samples.
4. Computes **mass depth** $m_i = \sum_k \rho_{db,k}\cdot\Delta z_k$
   (g/cm²) — cumulative mass per unit area down to each sample.
5. Saves a sibling CSV with the derived columns (DBD, Δz, mass depth) so you
   can audit the calculation.
6. Opens an interactive popup plot of **Excess Pb-210 vs mass depth**
   (log-x, mass depth on inverted y).
7. You click **two points** on the plot to bracket a regression interval.
   The script snaps each click to the nearest sample, selects every point
   between them, fits a linear regression to $\ln(\text{Excess Pb-210})$ vs
   mass depth, and computes:

   $$\mathrm{MAR} = -\frac{\lambda}{\text{slope}} \quad [\mathrm{g/(cm^2\cdot yr)}]$$

   where $\lambda = 0.03114\ \mathrm{yr^{-1}}$ (Pb-210 decay constant).
8. Repeats for as many intervals as you want; previous regression lines are
   shown faded on subsequent plots for context.
9. Saves a results CSV with the MAR, slope, R², and standard errors.

## Requirements

- Python 3.9+
- `pandas`, `numpy`, `matplotlib`
- `tkinter` (standard on Windows/macOS; `sudo apt install python3-tk` on Linux)

```powershell
pip install pandas numpy matplotlib
```

## Usage

```powershell
# Open file pickers for both CSVs
python mass_accumulation_rate.py

# Pass paths directly
python mass_accumulation_rate.py AgeModel.csv weights.csv

# Override the core label shown in plot titles
python mass_accumulation_rate.py AgeModel.csv weights.csv --core NBP2002KC11
```

If only one path is given, a dialog opens for the missing one.

## Input file requirements

### AgeModel CSV

(Produced by `210Pb_AgeModelPlot.ipynb`.) Required columns:

- `Center point of interval`
- `Top of interval (cm)`
- `Base of interval (cm)`
- `Excess Pb-210 (Bq/g)`
- `Pb-210 activity Uncertainty (Bq-g)`
- `calendar years pre year of core`

### Weights CSV

Required columns:

- `Center point of interval`
- `wet sediment weight (g)` — bulk wet mass of the volumetric sample
- `dry sediment weight (g)` — bulk dry mass after drying to constant weight

The two files are merged on `Center point of interval`. Rows that don't match
in both are dropped (with a warning).

## Output

Two CSVs are written next to the AgeModel CSV:

### `<AgeModel>_MassDepth_<YYYYMMDD>.csv` (derived table)

All the original AgeModel columns plus:

| column | meaning |
| --- | --- |
| `water content (wt fraction)` | $(M_w - M_d)/M_w$ |
| `DBD (g/cm3)` | dry bulk density (indirect method) |
| `effective thickness (cm)` | midpoint-to-midpoint Δz |
| `mass per area (g/cm2)` | `DBD × Δz` for that interval |
| `mass depth at base (g/cm2)` | cumulative mass per area down to the base of the interval |
| `mass depth at center (g/cm2)` | mass depth used as the plotting / regression x-coordinate (one-half of the increment back from the base) |

### `<AgeModel>_MAR_<YYYYMMDD>.csv` (results)

One row per clicked pair:

| column | meaning |
| --- | --- |
| `pair_index` | 1-based order you computed them in |
| `shallow_depth_cm` / `deep_depth_cm` | core depths bracketing the regression |
| `m_shallow` / `m_deep` | mass depths at those bounds (g/cm²) |
| `n_points` | how many samples went into the regression |
| `slope` ± `slope_se` | slope of $\ln(\text{Excess Pb-210})$ vs mass depth |
| `intercept` | regression intercept |
| `r_squared` | R² of the fit |
| `MAR_g_per_cm2_per_yr` | the MAR |
| `MAR_se_g_per_cm2_per_yr` | propagated standard error: $\sigma_{MAR} = \frac{\lambda}{\text{slope}^2}\,\sigma_{\text{slope}}$ |

## Step-by-step example

```powershell
python mass_accumulation_rate.py AgeModel_NBP2002KC11_20260614.csv NBP2002_KC11_weights.csv --core NBP2002KC11
```

1. The script prints the path to the derived `MassDepth` CSV so you can
   audit DBD and mass-depth values before clicking anything.
2. A window pops up titled `NBP2002KC11 - Click TWO points to bracket the
   regression interval (pair #1)`.
3. Click roughly at the shallowest sample you want in the regression, then at
   the deepest. Only the y (mass depth) coordinate is used; clicks snap to
   the nearest sample.
4. The window updates to show the red fitted line, the included samples
   highlighted, and the MAR in the title / legend. Close the window.
5. Terminal asks `Calculate another MAR interval? (yes/no):`. Type `yes` for
   another segment (e.g., to fit a shallower vs deeper population
   separately) or `no` to stop.
6. The MAR results CSV is written next to the AgeModel CSV.

## Method notes

- The interval-thickness convention used here treats each sample as
  representative of the depth range halfway to its neighbors. For a core
  with continuous 3 cm sampling, every Δz = 3 cm and this collapses to the
  naive method. For a gappy core, the unsampled stretches between samples
  are assigned to the nearest sample, so the mass-depth profile is
  continuous through the core. This is the convention assumed by Boldt
  et al. (2013) for 210Pb MAR work.
- The regression is fit to $\ln(\text{Excess Pb-210})$ vs mass depth.
  Samples with non-positive Excess Pb-210 are excluded from the fit
  (the log is undefined). They are still plotted on the profile (well, the
  positive ones are — non-positive are dropped from the log-x plot too).
- The MAR uncertainty propagates only the regression slope error
  ($\sigma_{MAR} = (\lambda/\text{slope}^2)\,\sigma_{\text{slope}}$).
  Activity-measurement uncertainty is shown as horizontal error bars on the
  plot but is not folded into the MAR standard error.
- The script assumes the **saturated** form of porosity (pore space = water,
  no air). If your samples had air voids or were not fully water-saturated,
  the indirect DBD will be biased. Use a direct (`dry mass / container
  volume`) calculation in that case.

## Troubleshooting

- **No window appears.** matplotlib can't open an interactive backend —
  usually a missing `tkinter`. On Linux: `sudo apt install python3-tk`.
- **`KeyError: ... missing required columns`.** Check the column names in
  both CSVs against the lists above. They must match exactly (including
  capitalization, spaces, and the `(g)` / `(cm)` suffixes).
- **"... rows in the AgeModel did not match a weights row".** Some
  `Center point of interval` values don't match between the two files
  (e.g., 1.5 in one and 1.50 in the other usually still matches because
  pandas treats them as equal floats; but if one has 1.5 and the other 2.0
  they won't). Open both CSVs and reconcile.
- **"Slope is non-negative".** The samples in your selected range don't
  show the expected decay of Excess Pb-210 with depth — either you picked
  too few points, picked across a discontinuity, or you've passed the
  background-corrected limit. Try a tighter or different range.
- **The deep samples land at much smaller mass depths than expected.** You
  may be looking at the legacy (option a) convention. This script always
  uses option b — the discrepancy is because option a was wrong, not because
  this is wrong.
