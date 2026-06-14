"""
Mass Accumulation Rate (MAR) calculator
---------------------------------------
Companion to 210Pb_AgeModelPlot.ipynb.

Reads an AgeModel CSV and a weights CSV, computes dry bulk density (indirect
Sanchez-Cabeza & Ruiz-Fernandez 2012 method) and mass depth for each interval,
then lets you click two points on a popup plot of Excess Pb-210 vs mass depth
to define a regression interval. MAR is computed as

    MAR = -lambda / slope_of_ln(Excess_Pb210)_vs_mass_depth      [g/(cm^2 yr)]

where lambda = 0.03114 yr^-1 (Pb-210 decay constant). Repeat for as many
intervals as you want; the results table is saved as a CSV.

Usage:
    python mass_accumulation_rate.py [AgeModel.csv] [weights.csv]

If either path is omitted, a file-picker dialog opens for it.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PB210_DECAY_CONSTANT = 0.03114  # yr^-1
RHO_WATER = 1.02   # g/cm^3 (seawater)
RHO_SOLID = 2.65   # g/cm^3 (quartz)

AGE_REQUIRED = [
    "Center point of interval",
    "Top of interval (cm)",
    "Base of interval (cm)",
    "Excess Pb-210 (Bq/g)",
    "Pb-210 activity Uncertainty (Bq-g)",
    "calendar years pre year of core",
]
WEIGHTS_REQUIRED = [
    "Center point of interval",
    "wet sediment weight (g)",
    "dry sediment weight (g)",
]


# ---------------------------------------------------------------------------
# File loading
# ---------------------------------------------------------------------------
def pick_csv_via_dialog(title: str) -> Path | None:
    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError:
        return None
    root = tk.Tk()
    root.withdraw()
    p = filedialog.askopenfilename(
        title=title,
        filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
    )
    root.destroy()
    return Path(p) if p else None


def load_csv(path: Path, required: list[str], label: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"{label} CSV is missing required columns: {missing}\nFile: {path}")
    return df


# ---------------------------------------------------------------------------
# Dry bulk density (indirect) + mass depth (option b: midpoint-to-midpoint)
# ---------------------------------------------------------------------------
def dry_bulk_density(Mw: pd.Series, Md: pd.Series,
                     rho_w: float = RHO_WATER, rho_s: float = RHO_SOLID) -> pd.Series:
    """Indirect dry bulk density from wet/dry masses.

    DBD = Md / (V_water + V_solid)
        = Md / ((Mw - Md)/rho_w + Md/rho_s)
    """
    return Md / ((Mw - Md) / rho_w + Md / rho_s)


def interval_thickness_midpoint(centers: np.ndarray,
                                tops: np.ndarray,
                                bases: np.ndarray) -> np.ndarray:
    """Option (b): each sample's effective thickness runs from the midpoint
    with its shallower neighbor to the midpoint with its deeper neighbor.
    Endpoints are clamped to the top of the first interval and the base of
    the last interval respectively. Assumes inputs are sorted shallow->deep.
    """
    n = len(centers)
    upper = np.empty(n)  # shallower bound (cm)
    lower = np.empty(n)  # deeper bound (cm)

    upper[0] = tops[0]
    lower[-1] = bases[-1]
    for i in range(1, n):
        midpt = 0.5 * (centers[i - 1] + centers[i])
        upper[i] = midpt
        lower[i - 1] = midpt

    dz = lower - upper
    return dz


def compute_mass_depth(df: pd.DataFrame) -> pd.DataFrame:
    """Add DBD, effective thickness, and cumulative mass-depth columns.

    Mass depth is the cumulative mass per unit area down to the BOTTOM of the
    sample's effective interval (g/cm^2). We also report the mass-depth at the
    center of the effective interval, which is the value used for plotting.
    """
    df = df.sort_values("Center point of interval").reset_index(drop=True)

    Mw = df["wet sediment weight (g)"].astype(float)
    Md = df["dry sediment weight (g)"].astype(float)
    df["water content (wt fraction)"] = (Mw - Md) / Mw
    df["DBD (g/cm3)"] = dry_bulk_density(Mw, Md)

    centers = df["Center point of interval"].to_numpy(float)
    tops = df["Top of interval (cm)"].to_numpy(float)
    bases = df["Base of interval (cm)"].to_numpy(float)
    dz = interval_thickness_midpoint(centers, tops, bases)
    df["effective thickness (cm)"] = dz

    incr = df["DBD (g/cm3)"].to_numpy() * dz
    df["mass per area (g/cm2)"] = incr
    df["mass depth at base (g/cm2)"] = np.cumsum(incr)
    df["mass depth at center (g/cm2)"] = df["mass depth at base (g/cm2)"] - 0.5 * incr
    return df


# ---------------------------------------------------------------------------
# Plotting + interactive MAR
# ---------------------------------------------------------------------------
def draw_profile(ax, df: pd.DataFrame, core_name: str, pair_index: int,
                 previous: list[dict]) -> None:
    valid = df["Excess Pb-210 (Bq/g)"] > 0  # log plot requires positives
    ax.errorbar(
        df.loc[valid, "Excess Pb-210 (Bq/g)"],
        df.loc[valid, "mass depth at center (g/cm2)"],
        xerr=df.loc[valid, "Pb-210 activity Uncertainty (Bq-g)"],
        fmt="o-", color="black", linewidth=1, markersize=5,
        capsize=3, ecolor="lightgrey", label="Excess Pb-210",
    )

    for prev in previous:
        m = np.linspace(prev["m_shallow"], prev["m_deep"], 50)
        a = np.exp(prev["intercept"] + prev["slope"] * m)
        ax.plot(a, m, "--", color="steelblue", alpha=0.6, linewidth=1)

    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Excess Pb-210 activity (Bq/g)", fontsize=12)
    ax.set_ylabel("Mass depth (g/cm$^2$)", fontsize=12)
    ax.set_title(
        f"{core_name} - Click TWO points to bracket the regression interval (pair #{pair_index})\n"
        f"close the window when finished",
        fontsize=11,
    )
    ax.grid(True, which="both", linestyle="-", linewidth=0.5, color="lightgray")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=9)


def snap_to_nearest(df: pd.DataFrame, mass_depth_click: float) -> int:
    """Return the index of the row with mass-depth closest to the click."""
    md = df["mass depth at center (g/cm2)"]
    return int((md - mass_depth_click).abs().idxmin())


def compute_mar_for_pair(df: pd.DataFrame, core_name: str, pair_index: int,
                        previous: list[dict]) -> dict | None:
    fig, ax = plt.subplots(figsize=(6, 9))
    draw_profile(ax, df, core_name, pair_index, previous)
    fig.tight_layout()

    print(f"\nPair #{pair_index}: click TWO points to bracket the regression interval "
          f"(only the y / mass-depth coordinate matters).")
    pts = plt.ginput(2, timeout=0, show_clicks=True)
    if len(pts) < 2:
        print("Did not receive two clicks; skipping.")
        plt.close(fig)
        return None

    i1 = snap_to_nearest(df, pts[0][1])
    i2 = snap_to_nearest(df, pts[1][1])
    lo, hi = sorted([i1, i2])
    sel = df.iloc[lo:hi + 1].copy()

    # Require positive activities for log regression and at least 2 points
    sel = sel[(sel["Excess Pb-210 (Bq/g)"] > 0)
              & sel["Excess Pb-210 (Bq/g)"].notna()
              & sel["mass depth at center (g/cm2)"].notna()]
    if len(sel) < 2:
        print(f"Need at least 2 points with positive Excess Pb-210; got {len(sel)}. Skipping.")
        plt.close(fig)
        return None

    m = sel["mass depth at center (g/cm2)"].to_numpy(float)
    a = sel["Excess Pb-210 (Bq/g)"].to_numpy(float)
    lna = np.log(a)

    # Linear regression ln(A) = intercept + slope * m
    slope, intercept = np.polyfit(m, lna, 1)
    yhat = intercept + slope * m
    ss_res = float(np.sum((lna - yhat) ** 2))
    ss_tot = float(np.sum((lna - lna.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    n = len(m)
    # Std error of slope (OLS)
    if n > 2 and np.var(m) > 0:
        sigma2 = ss_res / (n - 2)
        slope_se = float(np.sqrt(sigma2 / np.sum((m - m.mean()) ** 2)))
    else:
        slope_se = float("nan")

    if slope >= 0:
        print(f"Slope is non-negative ({slope:.4g}); cannot compute MAR for this interval.")
        plt.close(fig)
        return None

    mar = -PB210_DECAY_CONSTANT / slope
    mar_se = (PB210_DECAY_CONSTANT / slope ** 2) * slope_se if np.isfinite(slope_se) else float("nan")

    # Overlay the fitted line on the plot
    m_line = np.linspace(m.min(), m.max(), 50)
    a_line = np.exp(intercept + slope * m_line)
    ax.plot(a_line, m_line, "r-", linewidth=2.5,
            label=f"Pair #{pair_index}: MAR = {mar:.4g} g/cm$^2$/yr (R$^2$={r2:.3f})")
    ax.scatter(a, m, color="red", zorder=5, s=50, edgecolor="black")
    ax.set_title(f"{core_name} - pair #{pair_index}: MAR = {mar:.4g} g/cm$^2$/yr\n"
                 f"close this window to continue", fontsize=11)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), fontsize=9)
    fig.canvas.draw()

    print(f"  selected rows {lo}..{hi}  (n = {n})")
    print(f"  shallow depth: {sel['Center point of interval'].iloc[0]:.2f} cm  "
          f"(mass depth {m[0]:.3f} g/cm^2)")
    print(f"  deep    depth: {sel['Center point of interval'].iloc[-1]:.2f} cm  "
          f"(mass depth {m[-1]:.3f} g/cm^2)")
    print(f"  slope = {slope:.5g} +/- {slope_se:.2g}   R^2 = {r2:.4f}")
    print(f"  >>> MAR = {mar:.5g} +/- {mar_se:.2g}  g/(cm^2 yr)")

    plt.show()

    return {
        "pair_index": pair_index,
        "shallow_depth_cm": float(sel["Center point of interval"].iloc[0]),
        "deep_depth_cm":    float(sel["Center point of interval"].iloc[-1]),
        "m_shallow": float(m[0]),
        "m_deep":    float(m[-1]),
        "n_points":  int(n),
        "slope":     float(slope),
        "slope_se":  float(slope_se),
        "intercept": float(intercept),
        "r_squared": float(r2),
        "MAR_g_per_cm2_per_yr": float(mar),
        "MAR_se_g_per_cm2_per_yr": float(mar_se),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    parser = argparse.ArgumentParser(
        description="Interactive MAR calculator (Boldt et al. 2013 style)."
    )
    parser.add_argument("age_model", nargs="?", type=Path,
                        help="Path to AgeModel CSV. If omitted, a file dialog opens.")
    parser.add_argument("weights", nargs="?", type=Path,
                        help="Path to weights CSV. If omitted, a file dialog opens.")
    parser.add_argument("--core", default=None,
                        help="Core label for the plot title. Defaults to the age-model CSV stem.")
    args = parser.parse_args()

    age_path = args.age_model or pick_csv_via_dialog("Select AgeModel CSV")
    if age_path is None:
        print("No AgeModel CSV selected. Exiting."); return 1
    weights_path = args.weights or pick_csv_via_dialog("Select weights CSV")
    if weights_path is None:
        print("No weights CSV selected. Exiting."); return 1

    age_path = Path(age_path).expanduser().resolve()
    weights_path = Path(weights_path).expanduser().resolve()
    if not age_path.is_file():
        print(f"ERROR: file not found: {age_path}"); return 2
    if not weights_path.is_file():
        print(f"ERROR: file not found: {weights_path}"); return 2

    try:
        matplotlib.use("TkAgg")
    except Exception:
        pass

    age_df = load_csv(age_path, AGE_REQUIRED, "AgeModel")
    w_df = load_csv(weights_path, WEIGHTS_REQUIRED, "weights")

    merged = age_df.merge(
        w_df[["Center point of interval",
              "wet sediment weight (g)", "dry sediment weight (g)"]],
        on="Center point of interval", how="left",
    )
    if merged[["wet sediment weight (g)", "dry sediment weight (g)"]].isna().any().any():
        n_bad = merged[["wet sediment weight (g)", "dry sediment weight (g)"]].isna().any(axis=1).sum()
        print(f"WARNING: {n_bad} rows in the AgeModel did not match a weights row "
              f"on 'Center point of interval' and will be dropped.")
        merged = merged.dropna(subset=["wet sediment weight (g)", "dry sediment weight (g)"]).copy()

    df = compute_mass_depth(merged)

    # Save the derived table once so the user can audit DBD / mass depth
    date_tag = datetime.today().strftime("%Y%m%d")
    derived_path = age_path.with_name(f"{age_path.stem}_MassDepth_{date_tag}.csv")
    df.to_csv(derived_path, index=False)
    print(f"Derived DBD + mass-depth table saved -> {derived_path}")

    core_name = args.core or age_path.stem
    results: list[dict] = []
    pair_n = 1
    while True:
        result = compute_mar_for_pair(df, core_name, pair_n, results)
        if result is not None:
            results.append(result)
            pair_n += 1
        again = input("\nCalculate another MAR interval? (yes/no): ").strip().lower()
        if again not in ("y", "yes"):
            break

    if results:
        out = pd.DataFrame(results)
        out_path = age_path.with_name(f"{age_path.stem}_MAR_{date_tag}.csv")
        out.to_csv(out_path, index=False)
        print(f"\nSaved {len(out)} MAR interval(s) -> {out_path}")
        print(out.to_string(index=False))
    else:
        print("\nNo MAR intervals were calculated. Nothing saved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
