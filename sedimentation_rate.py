"""
Sedimentation rate calculator
-----------------------------
Standalone companion to 210Pb_AgeModelPlot.ipynb.

Usage:
    python sedimentation_rate.py [path/to/AgeModel.csv]

If no path is given, a file-picker dialog opens.

Workflow:
    1. The script opens an interactive plot window showing the age model
       (Excess Pb-210 vs depth, plus background activity).
    2. Click TWO points on the plot to set the DEPTH RANGE for a segment.
       Only the y-coordinate (depth) is used; each click is snapped to the
       nearest interval with a valid calendar year.
    3. Every dated point that falls within that range (not just the two you
       clicked) is fit with an ordinary least-squares regression of depth
       vs. calendar year. The sedimentation rate (cm/yr) is the (sign-
       corrected) slope of that fit, along with its standard error and R^2.
    4. A second window shows the fit itself (year vs depth, points + line)
       so you can check the regression before accepting it.
    5. Close both windows; the script asks whether to compute another
       segment.
    6. When you stop, every segment is written to
           <input_dir>/<input_basename>_SedimentationRates_<YYYYMMDD>.csv
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REQUIRED_COLS = [
    "Center point of interval",
    "calendar years pre year of core",
    "Excess Pb-210 (Bq/g)",
    "Top of interval (cm)",
    "Pb-210 activity Uncertainty (Bq-g)",
    "Averaged supported activity of Bi-214 and Pb-214 (Bq/g)",
    "Background activity uncertainty (Bq/g)",
]


def pick_csv_via_dialog() -> Path | None:
    """Open a native file-picker for the age-model CSV."""
    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError:
        return None
    root = tk.Tk()
    root.withdraw()
    path = filedialog.askopenfilename(
        title="Select AgeModel CSV",
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )
    root.destroy()
    return Path(path) if path else None


def load_age_model(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise KeyError(
            f"AgeModel CSV is missing required columns: {missing}\n"
            f"File: {csv_path}"
        )
    return df


def draw_age_model(ax, data: pd.DataFrame, core_name: str, seg_index: int,
                   previous_segments: list[dict]) -> None:
    valid_excess = ~data["Excess Pb-210 (Bq/g)"].isna()
    ax.plot(
        data.loc[valid_excess, "Excess Pb-210 (Bq/g)"],
        data.loc[valid_excess, "Center point of interval"],
        color="black", linewidth=1, zorder=2, label="Excess Activity",
    )

    yerr = np.abs(data["Center point of interval"] - data["Top of interval (cm)"])
    xerr = data["Pb-210 activity Uncertainty (Bq-g)"]
    for i in range(len(data)):
        x = data["Excess Pb-210 (Bq/g)"].iloc[i]
        y = data["Center point of interval"].iloc[i]
        if pd.isna(x):
            continue
        rect = patches.Rectangle(
            (x - xerr.iloc[i], y - yerr.iloc[i]),
            xerr.iloc[i] * 2, yerr.iloc[i] * 2,
            linewidth=0.5, edgecolor="grey", facecolor="lightgrey",
            alpha=0.5, zorder=1,
        )
        ax.add_patch(rect)

    valid_bg = ~data["Averaged supported activity of Bi-214 and Pb-214 (Bq/g)"].isna()
    ax.errorbar(
        data.loc[valid_bg, "Averaged supported activity of Bi-214 and Pb-214 (Bq/g)"],
        data.loc[valid_bg, "Center point of interval"],
        xerr=data.loc[valid_bg, "Background activity uncertainty (Bq/g)"],
        fmt="-", color="grey", label="Background Activity",
        capsize=5, linewidth=1, ecolor="darkgrey",
    )

    # Faded overlay of previously fit segments for context
    for prev in previous_segments:
        ax.plot(
            [prev["shallow_excess_pb210"], prev["deep_excess_pb210"]],
            [prev["shallow_depth_cm"],    prev["deep_depth_cm"]],
            "o--", color="steelblue", markersize=6, linewidth=1, alpha=0.6,
        )

    ax.set_xscale("log")
    ax.set_xlim(0.01, 10)
    ax.invert_yaxis()
    ax.set_xlabel("Bq/g", fontsize=12)
    ax.set_ylabel("Depth (cm)", fontsize=12)
    ax.set_title(
        f"{core_name} - Click TWO points to set the range (segment #{seg_index})\n"
        f"close the window when finished",
        fontsize=12,
    )
    ax.grid(True, which="both", linestyle="-", linewidth=0.5, color="lightgray")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=9)


def snap_to_nearest(lookup: pd.DataFrame, depth_click: float) -> tuple[float, float, float]:
    """Snap a clicked depth to the nearest interval with a valid age."""
    idx = (lookup["Center point of interval"] - depth_click).abs().idxmin()
    row = lookup.loc[idx]
    excess = row["Excess Pb-210 (Bq/g)"]
    return (
        float(row["Center point of interval"]),
        float(row["calendar years pre year of core"]),
        float(excess) if not pd.isna(excess) else float("nan"),
    )


def get_points_in_range(lookup: pd.DataFrame, d_min: float, d_max: float) -> pd.DataFrame:
    """All dated rows with depth in [d_min, d_max], sorted shallow -> deep."""
    mask = (lookup["Center point of interval"] >= d_min) & \
           (lookup["Center point of interval"] <= d_max)
    return lookup.loc[mask].sort_values("Center point of interval").reset_index(drop=True)


def fit_linear_regression(x: np.ndarray, y: np.ndarray) -> dict:
    """
    Ordinary least-squares fit of y as a function of x (numpy only, no
    scipy dependency). Returns slope, intercept, R^2, the standard error
    of the slope, and n. Standard error is NaN when n <= 2 (no residual
    degrees of freedom) or when x has zero variance.
    """
    n = len(x)
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    residuals = y - y_pred
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    sxx = np.sum((x - np.mean(x)) ** 2)
    if n > 2 and sxx > 0:
        mse = ss_res / (n - 2)
        std_err = np.sqrt(mse / sxx)
    else:
        std_err = float("nan")

    return {"slope": slope, "intercept": intercept, "r_squared": r_squared,
            "std_err": std_err, "n": n}


def _fmt_rate(rate: float, std_err: float) -> str:
    if np.isnan(std_err):
        return f"{rate:.4f} cm/yr"
    return f"{rate:.4f} \u00b1 {std_err:.4f} cm/yr"


def calculate_one_segment(data: pd.DataFrame, lookup: pd.DataFrame, core_name: str,
                          seg_index: int, previous_segments: list[dict]) -> dict | None:
    fig, ax = plt.subplots(figsize=(6, 9))
    draw_age_model(ax, data, core_name, seg_index, previous_segments)
    fig.tight_layout()

    print(f"\nSegment #{seg_index}: click TWO points to set the depth range "
          f"for the regression (only depth matters; x is ignored).")
    pts = plt.ginput(2, timeout=0, show_clicks=True)

    if len(pts) < 2:
        print("Did not receive two clicks; skipping.")
        plt.close(fig)
        return None

    d1, _, _ = snap_to_nearest(lookup, pts[0][1])
    d2, _, _ = snap_to_nearest(lookup, pts[1][1])
    if d1 > d2:
        d1, d2 = d2, d1

    subset = get_points_in_range(lookup, d1, d2)
    if len(subset) < 2:
        print(f"Only {len(subset)} dated point(s) fall between {d1:.2f} and "
              f"{d2:.2f} cm; need at least 2 to fit a line. Skipping.")
        plt.close(fig)
        return None
    if len(np.unique(subset["calendar years pre year of core"])) < 2:
        print("All points in range share the same calendar year; "
              "cannot fit a slope. Skipping.")
        plt.close(fig)
        return None
    if len(subset) == 2:
        print("Only 2 dated points in range: this reduces to the old "
              "point-to-point calculation, with no uncertainty estimate.")

    years = subset["calendar years pre year of core"].to_numpy(dtype=float)
    depths = subset["Center point of interval"].to_numpy(dtype=float)

    fit = fit_linear_regression(years, depths)
    # Depth is regressed on calendar year, but calendar year decreases with
    # depth (older = smaller year), so the raw slope is negative. Flip the
    # sign so a normal (deeper = older) profile reports a positive rate,
    # matching the old delta_depth / delta_years convention.
    rate = -fit["slope"]
    std_err_rate = fit["std_err"]  # magnitude is unaffected by the sign flip

    # Highlight every point actually used in the regression
    ax.plot(
        subset["Excess Pb-210 (Bq/g)"], subset["Center point of interval"],
        "ro", markersize=9, zorder=3,
        label=f"Segment #{seg_index}: n={fit['n']}, {rate:.4f} cm/yr",
    )
    ax.set_title(
        f"{core_name} - segment #{seg_index}: {_fmt_rate(rate, std_err_rate)} "
        f"(R\u00b2={fit['r_squared']:.3f}, n={fit['n']})\n"
        f"close both windows to continue",
        fontsize=11,
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=9)
    fig.canvas.draw()

    # Diagnostic plot: calendar year vs depth, with the fitted line, so the
    # fit quality is visible before you accept the segment.
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    ax2.scatter(years, depths, color="black", zorder=3, label="Dated points")
    year_line = np.linspace(years.min(), years.max(), 50)
    depth_line = fit["slope"] * year_line + fit["intercept"]
    ax2.plot(year_line, depth_line, color="red", linewidth=1.5, label="OLS fit")
    ax2.invert_yaxis()
    ax2.set_xlabel("Calendar year", fontsize=10)
    ax2.set_ylabel("Depth (cm)", fontsize=10)
    ax2.set_title(
        f"Segment #{seg_index} regression\n"
        f"rate = {_fmt_rate(rate, std_err_rate)}, R\u00b2 = {fit['r_squared']:.3f}",
        fontsize=10,
    )
    ax2.legend(fontsize=8)
    fig2.tight_layout()

    print(f"  depth range: {d1:.2f}-{d2:.2f} cm, {fit['n']} dated point(s) used")
    print(f"  >>> Sedimentation rate = {_fmt_rate(rate, std_err_rate)} "
          f"(R\u00b2 = {fit['r_squared']:.3f})")

    plt.show()  # blocks until both windows are closed

    return {
        "segment_index": seg_index,
        "shallow_depth_cm": d1,
        "deep_depth_cm": d2,
        "n_points": fit["n"],
        "sedimentation_rate_cm_per_yr": rate,
        "rate_std_err_cm_per_yr": std_err_rate,
        "r_squared": fit["r_squared"],
        "shallow_excess_pb210": subset["Excess Pb-210 (Bq/g)"].iloc[0],
        "deep_excess_pb210": subset["Excess Pb-210 (Bq/g)"].iloc[-1],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Interactive sedimentation rate calculator (popup-based, "
                     "regression over a clicked depth range)."
    )
    parser.add_argument(
        "csv", nargs="?", type=Path,
        help="Path to AgeModel CSV. If omitted, a file dialog opens."
    )
    parser.add_argument(
        "--core", default=None,
        help="Core label for the plot title. Defaults to the CSV filename stem."
    )
    args = parser.parse_args()

    csv_path = args.csv
    if csv_path is None:
        csv_path = pick_csv_via_dialog()
        if csv_path is None:
            print("No file selected. Exiting.")
            return 1

    csv_path = Path(csv_path).expanduser().resolve()
    if not csv_path.is_file():
        print(f"ERROR: file not found: {csv_path}")
        return 2

    # Prefer an interactive backend so the popup actually appears
    try:
        matplotlib.use("TkAgg")
    except Exception:
        pass  # fall back to whatever matplotlib picks

    data = load_age_model(csv_path)
    lookup = (data.dropna(subset=["Center point of interval",
                                  "calendar years pre year of core"])
                  .sort_values("Center point of interval")
                  .reset_index(drop=True))
    if lookup.empty:
        print("ERROR: no rows have both a depth and a calendar year.")
        return 3

    core_name = args.core or csv_path.stem
    results: list[dict] = []
    seg_n = 1

    while True:
        result = calculate_one_segment(data, lookup, core_name, seg_n, results)
        if result is not None:
            results.append(result)
            seg_n += 1
        again = input("\nCalculate another sedimentation rate? (yes/no): ").strip().lower()
        if again not in ("y", "yes"):
            break

    if results:
        df = pd.DataFrame(results)
        date_tag = datetime.today().strftime("%Y%m%d")
        out_path = csv_path.with_name(f"{csv_path.stem}_SedimentationRates_{date_tag}.csv")
        df.to_csv(out_path, index=False)
        print(f"\nSaved {len(df)} segment(s) -> {out_path}")
        print(df.to_string(index=False))
    else:
        print("\nNo sedimentation rates were calculated. Nothing saved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
