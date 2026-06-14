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
    2. Click TWO points on the plot. Only the y-coordinate (depth) is used;
       each click is snapped to the nearest interval with a valid calendar year.
    3. The sedimentation rate (cm/yr) is computed and shown on the plot.
    4. Close the window; the script asks whether to compute another pair.
    5. When you stop, every pair is written to
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


def draw_age_model(ax, data: pd.DataFrame, core_name: str, pair_index: int,
                   previous_pairs: list[dict]) -> None:
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

    # Faded overlay of previously calculated pairs for context
    for prev in previous_pairs:
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
        f"{core_name} - Click TWO points (pair #{pair_index})\n"
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


def calculate_one_pair(data: pd.DataFrame, lookup: pd.DataFrame, core_name: str,
                       pair_index: int, previous_pairs: list[dict]) -> dict | None:
    fig, ax = plt.subplots(figsize=(6, 9))
    draw_age_model(ax, data, core_name, pair_index, previous_pairs)
    fig.tight_layout()

    print(f"\nPair #{pair_index}: click TWO points in the popup window "
          f"(only depth matters; x is ignored).")
    pts = plt.ginput(2, timeout=0, show_clicks=True)

    if len(pts) < 2:
        print("Did not receive two clicks; skipping.")
        plt.close(fig)
        return None

    d1, y1, x1 = snap_to_nearest(lookup, pts[0][1])
    d2, y2, x2 = snap_to_nearest(lookup, pts[1][1])

    # Order shallow -> deep
    if d1 > d2:
        d1, d2 = d2, d1
        y1, y2 = y2, y1
        x1, x2 = x2, x1

    delta_depth = d2 - d1
    delta_years = y1 - y2  # shallower year (more recent) - deeper year (older)
    if delta_depth == 0 or delta_years == 0:
        print(f"Selected points collapse to a single interval "
              f"(Δdepth={delta_depth}, Δyears={delta_years}); cannot compute a rate.")
        plt.close(fig)
        return None

    rate = delta_depth / delta_years

    # Mark the selected segment, then keep the window open until the user closes it
    ax.plot([x1, x2], [d1, d2], "ro-", markersize=10, linewidth=2,
            label=f"Pair #{pair_index}: {rate:.4f} cm/yr")
    ax.set_title(
        f"{core_name} - pair #{pair_index}: {rate:.4f} cm/yr\n"
        f"close this window to continue",
        fontsize=12,
    )
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=9)
    fig.canvas.draw()

    print(f"  shallow: depth = {d1:.2f} cm, year = {y1:.1f}")
    print(f"  deep:    depth = {d2:.2f} cm, year = {y2:.1f}")
    print(f"  Δdepth = {delta_depth:.2f} cm,  Δyears = {delta_years:.2f} yr")
    print(f"  >>> Sedimentation rate = {rate:.4f} cm/yr  ({1/rate:.2f} yr/cm)")

    plt.show()  # blocks until the user closes the window

    return {
        "pair_index": pair_index,
        "shallow_depth_cm": d1,
        "shallow_year": y1,
        "shallow_excess_pb210": x1,
        "deep_depth_cm": d2,
        "deep_year": y2,
        "deep_excess_pb210": x2,
        "delta_depth_cm": delta_depth,
        "delta_years": delta_years,
        "sedimentation_rate_cm_per_yr": rate,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Interactive sedimentation rate calculator (popup-based)."
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
    pair_n = 1

    while True:
        result = calculate_one_pair(data, lookup, core_name, pair_n, results)
        if result is not None:
            results.append(result)
            pair_n += 1
        again = input("\nCalculate another sedimentation rate? (yes/no): ").strip().lower()
        if again not in ("y", "yes"):
            break

    if results:
        df = pd.DataFrame(results)
        date_tag = datetime.today().strftime("%Y%m%d")
        out_path = csv_path.with_name(f"{csv_path.stem}_SedimentationRates_{date_tag}.csv")
        df.to_csv(out_path, index=False)
        print(f"\nSaved {len(df)} pair(s) -> {out_path}")
        print(df.to_string(index=False))
    else:
        print("\nNo sedimentation rates were calculated. Nothing saved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
