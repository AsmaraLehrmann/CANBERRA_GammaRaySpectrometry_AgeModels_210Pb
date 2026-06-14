# Sedimentation Rate Calculator

Standalone Python script that lets you click two points on a 210Pb age model
and computes the sedimentation rate (cm/yr) between them. Designed as a
companion to `210Pb_AgeModelPlot.ipynb` — it reads the `AgeModel.csv` that
notebook produces.

## What it does

1. Loads an existing `AgeModel.csv`.
2. Opens an interactive popup window showing the age model
   (Excess Pb-210 + error rectangles + background activity, log-x, depth on y).
3. You click **two points** on the plot. Only the y-coordinate (depth) is
   used — each click is snapped to the nearest interval that has a valid
   calendar year, so you can't accidentally pick a depth with missing data.
4. The script:
   - orders the pair shallow → deep
   - computes `Δdepth = deep_depth − shallow_depth` (cm)
   - computes `Δyears = shallow_year − deep_year` (yr)
   - reports `rate = Δdepth / Δyears` (cm/yr) and its inverse (yr/cm)
   - draws the selected segment in red on the plot
5. After you close the window, it asks whether to calculate another pair.
   Earlier pairs are shown in faded blue on subsequent plots for context.
6. When you finish, every pair is saved to a CSV next to the input file.

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
blank are skipped automatically when snapping clicks.

## Output

A CSV is written next to the input file:

```
<input_basename>_SedimentationRates_<YYYYMMDD>.csv
```

Columns:

| column | meaning |
| --- | --- |
| `pair_index` | 1-based index in the order you computed the pairs |
| `shallow_depth_cm` | depth of the shallower snapped point (cm) |
| `shallow_year` | calendar year at that depth |
| `shallow_excess_pb210` | Excess Pb-210 value used to draw the red point |
| `deep_depth_cm` | depth of the deeper snapped point (cm) |
| `deep_year` | calendar year at that depth |
| `deep_excess_pb210` | Excess Pb-210 value used to draw the red point |
| `delta_depth_cm` | `deep_depth_cm − shallow_depth_cm` |
| `delta_years` | `shallow_year − deep_year` |
| `sedimentation_rate_cm_per_yr` | `delta_depth_cm / delta_years` |

Nothing is written if you exit without computing any pairs.

## Step-by-step example

```powershell
python sedimentation_rate.py "C:\Users\aalehrma\Desktop\210Pb\runs\MB1901\AgeModel_MB1901_20260614.csv" --core MB1901
```

1. A window pops up titled `MB1901 - Click TWO points (pair #1)`.
2. Click roughly on the depths you want. The clicks snap to the nearest
   intervals — that's intentional and means you don't have to be precise.
3. The plot updates to show the red segment and the rate in the title /
   legend; close the window.
4. Terminal asks `Calculate another sedimentation rate? (yes/no):` — type
   `yes` to do another pair (the previous one stays drawn faintly for
   reference) or `no` to stop.
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
- **The pair collapses to a single interval.** Both clicks snapped to the
  same depth, or the two snapped depths happen to share the same calendar
  year. Pick points that are farther apart.

## Notes

- The script does not modify the input CSV.
- It writes a sibling CSV with a date tag; running it again the same day
  overwrites that file.
- The plot window is independent — you can pan/zoom with the matplotlib
  toolbar before clicking. Clicks during pan/zoom mode won't register
  (matplotlib default behavior); exit those modes before clicking the two
  points.
