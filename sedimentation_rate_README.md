# Sedimentation Rate Calculator

Standalone Python script that lets you click a depth range on a 210Pb age
model and fits a linear regression through every dated point in that range
to compute the sedimentation rate (cm/yr). Designed as a companion to
`210Pb_AgeModelPlot.ipynb` — it reads the `AgeModel.csv` that notebook
produces.

## What it does

1. Loads an existing `AgeModel.csv`.
2. Opens an interactive popup window showing the age model
   (Excess Pb-210 + error rectangles + background activity, log-x, depth on y).
3. You click **two points** on the plot. Only the y-coordinate (depth) is
   used — each click is snapped to the nearest interval that has a valid
   calendar year. The two snapped depths become the **boundaries of a
   range**, not the only data used.
4. The script:
   - collects every dated interval whose depth falls inside that range
   - fits an ordinary least-squares (OLS) regression of depth against
     calendar year across all of those points (see [Equations](#equations)
     below)
   - reports the rate, its standard error, and R² (goodness of fit)
   - highlights every point used (not just the two clicked) in red on the
     age-model plot
   - opens a **second window** showing calendar year vs. depth for just
     the points used, with the fitted line drawn through them, so you can
     see the fit quality before trusting the number
5. After you close both windows, it asks whether to calculate another
   segment. Earlier segments are shown in faded blue on subsequent plots
   for context.
6. When you finish, every segment is saved to a CSV next to the input file.

## Equations

The age model gives each depth interval a calendar year via the constant
initial concentration (CIC) ²¹⁰Pb model:

```
Age = (1 / λ) × ln(A_surface / A_depth)
```

where λ = 0.03114 yr⁻¹ (²¹⁰Pb decay constant, ln(2)/22.3) and A is excess
²¹⁰Pb activity. `calendar years pre year of core` in the CSV is the
calendar year derived from that age.

**Regression.** For the points inside your selected range, the script fits
depth as a linear function of calendar year:

```
depth = a + b × year
```

using ordinary least squares (`numpy.polyfit`, degree 1). `b` is the raw
slope in cm per calendar year.

**Sign correction.** Calendar year decreases with depth (older = smaller
year), so the raw slope `b` comes out negative for a normal profile. The
reported rate flips the sign so a normal (deeper = older) profile reads
positive, matching the old point-to-point convention:

```
sedimentation_rate_cm_per_yr = −b
```

**Goodness of fit (R²).**

```
R² = 1 − (Σ(depth_i − depth_predicted_i)²) / (Σ(depth_i − depth_mean)²)
```

**Standard error of the rate.** Only defined when there are more than 2
points in the range (residual degrees of freedom > 0):

```
MSE = Σ(depth_i − depth_predicted_i)² / (n − 2)
SE(b) = √(MSE / Σ(year_i − year_mean)²)
```

`SE(b)` and `SE(rate)` are the same magnitude — the sign flip above
doesn't change an error bar. With exactly 2 points in range, `SE` is
undefined (no residual degrees of freedom) and is reported as blank/NaN;
this is mathematically identical to the old two-point calculation, just
computed via the regression code path.

This mirrors the standard published form for a single-segment fit,
`S = λ × (−m)`, where `m` is the slope of depth regressed on ln(activity)
directly. The script here regresses depth against calendar year instead of
raw ln(activity), since the CSV already carries year values computed from
the CIC model. The two are equivalent when the CIC model holds throughout
the fitted range — the regression here is just working one algebraic step
downstream of the raw activity data.

## Requirements

- Python 3.9+
- Packages: `pandas`, `numpy`, `matplotlib`
- `tkinter` (bundled with standard Python on Windows / macOS; on Linux install
  via your package manager, e.g. `sudo apt install python3-tk`)

Install the Python packages if needed:

```powershell
pip install pandas numpy matplotlib
```

## Usage

From a terminal in the project folder:

```powershell
# Open a file-picker dialog to choose the CSV
python sedimentation_rate.py

# Pass the CSV path directly
python sedimentation_rate.py path\to\AgeModel.csv

# Pass the CSV and override the core label shown in the plot title
python sedimentation_rate.py path\to\AgeModel.csv --core MB1901
```

If `--core` is omitted, the CSV filename stem is used as the plot title label.

## Input file requirements

The script reads a CSV with these columns (the names match what
`210Pb_AgeModelPlot.ipynb` writes):

- `Center point of interval`
- `calendar years pre year of core`
- `Excess Pb-210 (Bq/g)`
- `Top of interval (cm)`
- `Pb-210 activity Uncertainty (Bq-g)`
- `Averaged supported activity of Bi-214 and Pb-214 (Bq/g)`
- `Background activity uncertainty (Bq/g)`

If any are missing, the script stops with a clear error listing which.
Rows where `Center point of interval` or `calendar years pre year of core` is
blank are skipped automatically when snapping clicks or building a range.

## Output

A CSV is written next to the input file:

```
<input_basename>_SedimentationRates_<YYYYMMDD>.csv
```

Columns:

| column | meaning |
| --- | --- |
| `segment_index` | 1-based index in the order you computed the segments |
| `shallow_depth_cm` | shallow boundary of the clicked range (cm) |
| `deep_depth_cm` | deep boundary of the clicked range (cm) |
| `n_points` | number of dated intervals used in the regression |
| `sedimentation_rate_cm_per_yr` | sign-corrected regression slope, `−b` |
| `rate_std_err_cm_per_yr` | standard error of the rate; blank/NaN if `n_points == 2` |
| `r_squared` | R² of the depth-vs-year fit |
| `shallow_excess_pb210` | Excess Pb-210 at the shallowest point used |
| `deep_excess_pb210` | Excess Pb-210 at the deepest point used |

Nothing is written if you exit without computing any segments.

## Step-by-step example

```powershell
python sedimentation_rate.py "C:\Users\aalehrma\Desktop\210Pb\runs\MB1901\AgeModel_MB1901_20260614.csv" --core MB1901
```

1. A window pops up titled `MB1901 - Click TWO points to set the range (segment #1)`.
2. Click roughly on the shallow and deep bounds of the interval you want.
   The clicks snap to the nearest dated intervals — that's intentional and
   means you don't have to be precise.
3. The plot updates to show every point used (not just your two clicks) in
   red, with the rate, R², and n in the title. A second, smaller window
   opens showing the year-vs-depth fit itself.
4. Close both windows. Terminal asks
   `Calculate another sedimentation rate? (yes/no):` — type `yes` to do
   another segment (the previous one stays drawn faintly for reference) or
   `no` to stop.
5. The summary table prints and a CSV is written:
   `AgeModel_MB1901_20260614_SedimentationRates_<today>.csv`.

## Troubleshooting

- **No window appears.** Matplotlib couldn't open an interactive backend.
  The script tries to force `TkAgg`, which needs `tkinter`. On Linux,
  install `python3-tk` and rerun. On macOS, use a Python build that
  ships with Tk (the python.org installer or `brew install python-tk`).
- **`KeyError: AgeModel CSV is missing required columns: [...]`.** The CSV
  doesn't have the columns listed under "Input file requirements" above.
  Either you selected the wrong CSV or you're using an older format.
  Re-run `210Pb_AgeModelPlot.ipynb` to regenerate the AgeModel CSV.
- **The clicks don't behave as expected.** The x-axis is ignored — only the
  vertical position of each click matters. Click anywhere at the correct
  depth.
- **"Did not receive two clicks; skipping."** You closed the window before
  clicking twice. Answer `yes` to the next prompt and try again.
- **"Only N dated point(s) fall between ... need at least 2."** Your range
  is too narrow — widen it so it spans at least two dated intervals.
- **"All points in range share the same calendar year."** The regression
  needs variation in calendar year to fit a slope; widen the range.
- **`rate_std_err_cm_per_yr` is blank.** You selected a range with exactly
  2 dated points, so there are no residual degrees of freedom to estimate
  uncertainty from. The rate itself is still valid — it's identical to the
  old point-to-point calculation — there's just no error bar attached.

## Notes

- The script does not modify the input CSV.
- It writes a sibling CSV with a date tag; running it again the same day
  overwrites that file.
- The plot window is independent — you can pan/zoom with the matplotlib
  toolbar before clicking. Clicks during pan/zoom mode won't register
  (matplotlib default behavior); exit those modes before clicking the two
  boundary points.
- Wider ranges give you more points and a more stable rate, but they also
  average over more of the profile — if there's a real slope break inside
  your range, the regression will blend both segments into one number
  rather than showing you where the break is. Pick ranges you have reason
  to believe are internally consistent (a facies unit, a segment already
  identified as linear on the raw plot), not just "as wide as possible."
