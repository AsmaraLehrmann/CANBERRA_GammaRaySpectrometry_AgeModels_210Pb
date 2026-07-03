"""
210Pb decay / MAR curve plotter
-------------------------------
Interactive plot customizer for 210Pb work. Companion to
210Pb_AgeModelPlot.ipynb, sedimentation_rate.py, mass_accumulation_rate.py.

Two modes:
    decay  - Excess Pb-210 (log x) vs depth (cm), inverted y
    mar    - Excess Pb-210 (log x) vs mass depth (g/cm^2), inverted y

Inputs:
    decay : AgeModel CSV (produced by the notebook)
    mar   : MassDepth CSV (produced by mass_accumulation_rate.py)

Interactive controls (in one matplotlib window):
    - Title text box: type to retitle the plot
    - X min / X max / Y top / Y bot text boxes: type to set axis limits
    - "Add year"      button: next click on the plot adds a calendar-year (CE)
                              annotation at the nearest sample
    - "Mark missing"  button: next click marks the nearest sample's interval
                              with a brown band (for undetectable radioisotopes)
    - "Undo year" / "Undo missing": pop the most recent annotation/marker
    - "Save PDF": writes a timestamped PDF next to the input CSV
    - "Done": close the window

Usage:
    python curve_plotter.py [--mode decay|mar] [CSV] [--core LABEL]

If --mode is omitted, the script prompts. If CSV is omitted, a file dialog opens.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, ScalarFormatter
from matplotlib.widgets import Button, TextBox
import numpy as np
import pandas as pd

REQUIRED_BASE = [
    "Center point of interval",
    "Top of interval (cm)",
    "Base of interval (cm)",
    "Excess Pb-210 (Bq/g)",
    "Pb-210 activity Uncertainty (Bq-g)",
    "calendar years pre year of core",
]
BACKGROUND_COLS = [
    "Averaged supported activity of Bi-214 and Pb-214 (Bq/g)",
    "Background activity uncertainty (Bq/g)",
]
MAR_Y_COL = "mass depth at center (g/cm2)"


def pick_csv(title: str) -> Path | None:
    try:
        import tkinter as tk
        from tkinter import filedialog
    except ImportError:
        return None
    root = tk.Tk()
    root.withdraw()
    p = filedialog.askopenfilename(
        title=title, filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
    )
    root.destroy()
    return Path(p) if p else None


class CurveEditor:
    def __init__(self, df: pd.DataFrame, mode: str, core_name: str, save_dir: Path):
        self.df = df
        self.mode = mode  # 'decay' or 'mar'
        self.core_name = core_name
        self.save_dir = Path(save_dir)
        self.year_labels: list = []
        self.missing_spans: list = []
        self.current_action: str | None = None  # None | 'year' | 'missing'

        self.fig = plt.figure(figsize=(3, 10))
        self.ax = self.fig.add_axes([0.15, 0.44, 0.78, 0.48])
        self._widget_axes: list = []
        self._draw_base()
        self._add_controls()

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)

    # ---------- plotting ----------
    def _ycol(self) -> str:
        return MAR_Y_COL if self.mode == "mar" else "Center point of interval"

    def _ylabel(self) -> str:
        return "Mass depth (g/cm$^2$)" if self.mode == "mar" else "Depth (cm)"

    def _y_errors(self, mask: pd.Series) -> np.ndarray:
        """2xN array of (lower, upper) y errors representing each sample's interval span."""
        df = self.df.loc[mask]
        if self.mode == "mar":
            # Each sample's mass-depth band ~ mass_per_area centered on mass_depth_at_center.
            if "mass per area (g/cm2)" in df.columns:
                half = df["mass per area (g/cm2)"].to_numpy(float) / 2.0
            else:
                spacing = self.df[self._ycol()].diff().abs().median()
                fallback = (spacing if spacing and not pd.isna(spacing) else 1.0) / 2.0
                half = np.full(len(df), fallback)
            return np.vstack([half, half])
        # decay mode: cm intervals (asymmetric, in case center is off-midpoint)
        center = df["Center point of interval"].to_numpy(float)
        top    = df["Top of interval (cm)"].to_numpy(float)
        base   = df["Base of interval (cm)"].to_numpy(float)
        lower = np.clip(center - top, 0, None)
        upper = np.clip(base - center, 0, None)
        return np.vstack([lower, upper])

    def _draw_base(self) -> None:
        df = self.df
        y = df[self._ycol()]
        valid = df["Excess Pb-210 (Bq/g)"] > 0
        yerr = self._y_errors(valid)  # shape (2, N): [lower, upper]
        xvals = df.loc[valid, "Excess Pb-210 (Bq/g)"].to_numpy(float)
        yvals = y[valid].to_numpy(float)
        xerr  = df.loc[valid, "Pb-210 activity Uncertainty (Bq-g)"].to_numpy(float)

        # Light-grey uncertainty boxes (x: activity error, y: sample interval span)
        for xi, yi, xe, ylo, yup in zip(xvals, yvals, xerr, yerr[0], yerr[1]):
            if not (np.isfinite(xi) and np.isfinite(yi) and xi - xe > 0):
                continue
            rect = patches.Rectangle(
                (xi - xe, yi - ylo),
                width=2 * xe,
                height=ylo + yup,
                linewidth=0.5,
                edgecolor="grey",
                facecolor="lightgrey",
                alpha=0.5,
                zorder=1,
            )
            self.ax.add_patch(rect)

        # Central line + markers (no whiskers)
        self.ax.plot(
            xvals, yvals, "-", color="black", linewidth=1,
            zorder=3, label="Excess Activity",
        )
        self.ax.plot(
            xvals, yvals, "o", color="black", markersize=4, zorder=4,
        )
        if all(c in df.columns for c in BACKGROUND_COLS):
            bg = df["Averaged supported activity of Bi-214 and Pb-214 (Bq/g)"]
            bgerr = df["Background activity uncertainty (Bq/g)"]
            m = bg.notna() & (bg > 0)
            if m.any():
                self.ax.errorbar(
                    bg[m], y[m], xerr=bgerr[m],
                    fmt="-", color="grey", capsize=3, linewidth=1,
                    ecolor="darkgrey", label="Background activity", zorder=2,
                )
        self.ax.set_xscale("log")
        self.ax.invert_yaxis()
        self.ax.set_xlabel("Excess Pb-210 (Bq/g)", fontsize=12)
        self.ax.set_ylabel(self._ylabel(), fontsize=12)
        # Default x-min at 0.01 Bq/g
        self.ax.set_xlim(left=0.01)
        # Only label powers of 10 on the log x-axis (avoids crowding on narrow figures)
        self.ax.xaxis.set_minor_formatter(NullFormatter())
        # Plain decimal labels (0.01, 0.1, 1) instead of 10^-2, 10^-1, ...
        major_fmt = ScalarFormatter()
        major_fmt.set_scientific(False)
        major_fmt.set_useOffset(False)
        self.ax.xaxis.set_major_formatter(major_fmt)
        kind = "MAR Curve" if self.mode == "mar" else "210Pb Decay Curve"
        self.title = f"{self.core_name} {kind}"
        self.ax.set_title(self.title, fontsize=14)
        self.ax.grid(True, which="both", linestyle="-", linewidth=0.5, color="lightgray")
        # Legend below the x-axis label, single row, out of the data area
        self.ax.legend(
            loc="upper center", bbox_to_anchor=(0.5, -0.13),
            ncol=2, fontsize=9, frameon=False,
        )

    # ---------- widgets ----------
    def _add_controls(self) -> None:
        # Title
        ax_title = self.fig.add_axes([0.20, 0.32, 0.70, 0.04])
        self._widget_axes.append(ax_title)
        self.tb_title = TextBox(ax_title, "Title ", initial=self.title)
        self.tb_title.on_submit(self._set_title)

        # Axis text boxes (with current values pre-filled)
        x0, x1 = self.ax.get_xlim()
        y0, y1 = self.ax.get_ylim()  # inverted: y0 > y1
        y_top, y_bot = min(y0, y1), max(y0, y1)
        ax_xmin = self.fig.add_axes([0.10, 0.255, 0.15, 0.04])
        ax_xmax = self.fig.add_axes([0.32, 0.255, 0.15, 0.04])
        ax_ytop = self.fig.add_axes([0.55, 0.255, 0.15, 0.04])
        ax_ybot = self.fig.add_axes([0.77, 0.255, 0.15, 0.04])
        self._widget_axes += [ax_xmin, ax_xmax, ax_ytop, ax_ybot]
        self.tb_xmin = TextBox(ax_xmin, "X min ", initial=f"{x0:g}")
        self.tb_xmax = TextBox(ax_xmax, "X max ", initial=f"{x1:g}")
        self.tb_ytop = TextBox(ax_ytop, "Y top ", initial=f"{y_top:g}")
        self.tb_ybot = TextBox(ax_ybot, "Y bot ", initial=f"{y_bot:g}")
        self.tb_xmin.on_submit(lambda _t: self._set_xlim())
        self.tb_xmax.on_submit(lambda _t: self._set_xlim())
        self.tb_ytop.on_submit(lambda _t: self._set_ylim())
        self.tb_ybot.on_submit(lambda _t: self._set_ylim())

        # Missing intervals text box (typed, comma-separated cm ranges)
        ax_missing = self.fig.add_axes([0.20, 0.19, 0.50, 0.04])
        self._widget_axes.append(ax_missing)
        self.tb_missing = TextBox(
            ax_missing, "Missing (cm) ",
            initial="",
        )
        self.tb_missing.on_submit(self._set_missing)

        # Figure-size text box, on the same row to the right of "Missing"
        w_in, h_in = self.fig.get_size_inches()
        ax_size = self.fig.add_axes([0.82, 0.19, 0.11, 0.04])
        self._widget_axes.append(ax_size)
        self.tb_size = TextBox(
            ax_size, "Fig WxH (in) ",
            initial=f"{w_in:g}x{h_in:g}",
        )
        self.tb_size.on_submit(lambda _t: self._set_size())

        # Buttons row
        bw, bh, by = 0.18, 0.05, 0.10
        cols = [0.06, 0.30, 0.54, 0.78]
        labels_actions = [
            ("Add year",  lambda e: self._set_action("year")),
            ("Undo year", lambda e: self._undo_year()),
            ("Save PDF",  lambda e: self._save()),
            ("Done",      lambda e: plt.close(self.fig)),
        ]
        self._buttons = []
        for col, (label, cb) in zip(cols, labels_actions):
            ax_b = self.fig.add_axes([col, by, bw, bh])
            self._widget_axes.append(ax_b)
            btn = Button(ax_b, label)
            btn.on_clicked(cb)
            self._buttons.append(btn)

        # Status bar
        ax_status = self.fig.add_axes([0.04, 0.02, 0.92, 0.05]); ax_status.axis("off")
        self._widget_axes.append(ax_status)
        self.status = ax_status.text(
            0.5, 0.5,
            "Ready. Examples: Missing (cm) = '1-2, 5-8'.  Add year -> click sample.",
            ha="center", va="center", fontsize=10, color="darkblue"
        )

    # ---------- callbacks ----------
    def _set_status(self, msg: str) -> None:
        self.status.set_text(msg); self.fig.canvas.draw_idle()

    def _set_action(self, action: str) -> None:
        self.current_action = action
        self._set_status("Click on the plot near a sample to add a calendar-year (CE) label.")

    def _set_title(self, text: str) -> None:
        self.title = text
        self.ax.set_title(text, fontsize=14)
        self.fig.canvas.draw_idle()

    def _set_xlim(self) -> None:
        try:
            x0 = float(self.tb_xmin.text); x1 = float(self.tb_xmax.text)
            if x0 <= 0 or x1 <= 0:
                raise ValueError("log x requires positive limits")
            if x0 == x1:
                raise ValueError("xmin == xmax")
            self.ax.set_xlim(x0, x1); self.fig.canvas.draw_idle()
            self._set_status(f"x-axis set to [{x0:g}, {x1:g}].")
        except Exception as e:
            self._set_status(f"Bad x limit: {e}")

    def _set_size(self) -> None:
        text = self.tb_size.text.strip().lower()
        # Accept "8x10", "8 x 10", "8,10", or "8 10"
        for sep in ("x", ",", " "):
            if sep in text:
                parts = [p for p in text.split(sep) if p.strip()]
                if len(parts) == 2:
                    try:
                        w = float(parts[0]); h = float(parts[1])
                        if w <= 0 or h <= 0:
                            raise ValueError("dimensions must be positive")
                        self.fig.set_size_inches(w, h, forward=True)
                        self.fig.canvas.draw_idle()
                        self._set_status(f"Figure resized to {w:g} x {h:g} in.")
                        return
                    except Exception as e:
                        self._set_status(f"Bad figure size: {e}")
                        return
        self._set_status("Bad figure size: use WxH (e.g. 8x10).")

    def _set_ylim(self) -> None:
        try:
            yt = float(self.tb_ytop.text); yb = float(self.tb_ybot.text)
            if yt == yb:
                raise ValueError("Y top == Y bot")
            # Invert: top (smaller value) appears at the top of the inverted axis
            self.ax.set_ylim(max(yt, yb), min(yt, yb)); self.fig.canvas.draw_idle()
            self._set_status(f"y-axis set to top={min(yt,yb):g}, bottom={max(yt,yb):g}.")
        except Exception as e:
            self._set_status(f"Bad y limit: {e}")

    def _snap_to_sample(self, y_click: float) -> pd.Series:
        sub = self.df.dropna(subset=[self._ycol()]).reset_index(drop=True)
        idx = (sub[self._ycol()] - y_click).abs().idxmin()
        return sub.loc[idx]

    def _on_click(self, event) -> None:
        if event.inaxes is not self.ax or self.current_action is None:
            return
        if event.ydata is None:
            return
        row = self._snap_to_sample(event.ydata)
        if self.current_action == "year":
            year = row.get("calendar years pre year of core", np.nan)
            if pd.isna(year):
                self._set_status("Nearest sample has no calendar year; not labeled.")
            else:
                x = row["Excess Pb-210 (Bq/g)"]
                y_val = row[self._ycol()]
                if pd.isna(x) or x <= 0:
                    x = self.ax.get_xlim()[1] * 0.7
                year_int = int(np.ceil(float(year)))
                txt = self.ax.text(x * 1.15, y_val, f"{year_int} CE",
                                    fontsize=11, color="black", va="center", zorder=4)
                self.year_labels.append(txt)
                self._set_status(f"Added '{year_int} CE' at y={y_val:.2f}.")
        self.current_action = None
        self.fig.canvas.draw_idle()

    # ---------- missing intervals (typed) ----------
    def _parse_missing(self, text: str) -> list[tuple[float, float]]:
        """Parse e.g. '1-2, 5-8 cm; 12-15' into [(1,2),(5,8),(12,15)]."""
        out: list[tuple[float, float]] = []
        if not text or not text.strip():
            return out
        cleaned = text.lower().replace("cm", "")
        # split on commas or semicolons
        parts = [p.strip() for p in cleaned.replace(";", ",").split(",") if p.strip()]
        for part in parts:
            # accept 1-2, 1 - 2, 1 to 2
            tok = part.replace("to", "-")
            if "-" not in tok:
                raise ValueError(f"can't parse '{part}' (expected 'lo-hi')")
            lo_s, hi_s = tok.split("-", 1)
            lo = float(lo_s.strip()); hi = float(hi_s.strip())
            if lo > hi:
                lo, hi = hi, lo
            out.append((lo, hi))
        return out

    def _cm_to_y(self, depth_cm: float) -> float:
        """Convert a core depth in cm to the plot's y-coordinate.

        In decay mode this is the identity. In MAR mode, interpolate using
        the (center_cm, mass_depth_at_center) curve with (0, 0) anchored at
        the surface.
        """
        if self.mode != "mar":
            return depth_cm
        df = self.df.dropna(subset=["Center point of interval", MAR_Y_COL]).copy()
        df = df.sort_values("Center point of interval")
        xs = np.concatenate([[0.0], df["Center point of interval"].to_numpy(float)])
        ys = np.concatenate([[0.0], df[MAR_Y_COL].to_numpy(float)])
        return float(np.interp(depth_cm, xs, ys))

    def _set_missing(self, text: str) -> None:
        # Clear any existing missing spans first
        for s in self.missing_spans:
            try:
                s.remove()
            except Exception:
                pass
        self.missing_spans = []

        try:
            intervals = self._parse_missing(text)
        except Exception as e:
            self._set_status(f"Bad missing list: {e}")
            self.fig.canvas.draw_idle()
            return

        for lo_cm, hi_cm in intervals:
            top_y = self._cm_to_y(lo_cm)
            bot_y = self._cm_to_y(hi_cm)
            span = self.ax.axhspan(top_y, bot_y, color="brown", alpha=0.35, zorder=0)
            self.missing_spans.append(span)

        if intervals:
            joined = ", ".join(f"{lo:g}-{hi:g}" for lo, hi in intervals)
            self._set_status(f"Marked missing intervals: {joined} cm.")
        else:
            self._set_status("Cleared all missing-interval markers.")
        self.fig.canvas.draw_idle()

    def _undo_year(self) -> None:
        if not self.year_labels:
            self._set_status("No year labels to undo."); return
        self.year_labels.pop().remove(); self.fig.canvas.draw_idle()
        self._set_status("Removed last year label.")

    def _save(self) -> None:
        # Make sure the legend reflects any added "missing" patch
        _handles, labels = self.ax.get_legend_handles_labels()
        if self.missing_spans and "No radioisotopes detected" not in labels:
            self.missing_spans[0].set_label("No radioisotopes detected")
            self.ax.legend(
                loc="upper center", bbox_to_anchor=(0.5, -0.13),
                ncol=2, fontsize=9, frameon=False,
            )

        date_tag = datetime.today().strftime("%Y%m%d")
        kind = "MARCurve" if self.mode == "mar" else "DecayCurve"
        fname = f"{self.core_name}_{kind}_{date_tag}.pdf"
        out = self.save_dir / fname

        # Hide every widget (text boxes, buttons, status) so the saved PDF
        # contains only the main plot and its legend/title/labels.
        hidden = [a for a in self._widget_axes if a.get_visible()]
        for a in hidden:
            a.set_visible(False)
        try:
            # bbox_inches='tight' will now use only the visible artists (the plot).
            self.fig.savefig(out, dpi=300, bbox_inches="tight")
        finally:
            for a in hidden:
                a.set_visible(True)
            self.fig.canvas.draw_idle()
        self._set_status(f"Saved -> {out}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Interactive 210Pb decay / MAR curve plotter.")
    parser.add_argument("--mode", choices=["decay", "mar"], default=None)
    parser.add_argument("csv", nargs="?", type=Path,
                        help="Path to the CSV (AgeModel for decay, MassDepth for mar).")
    parser.add_argument("--core", default=None, help="Core label for the plot title.")
    args = parser.parse_args()

    mode = args.mode
    if mode is None:
        mode = input("Curve type - 'decay' or 'mar'? ").strip().lower()
        if mode not in ("decay", "mar"):
            print("Mode must be 'decay' or 'mar'."); return 1

    csv_path = args.csv
    if csv_path is None:
        title = "Select AgeModel CSV" if mode == "decay" else "Select MassDepth CSV"
        csv_path = pick_csv(title)
    if csv_path is None:
        print("No CSV selected. Exiting."); return 1
    csv_path = Path(csv_path).expanduser().resolve()
    if not csv_path.is_file():
        print(f"ERROR: file not found: {csv_path}"); return 2

    try:
        matplotlib.use("TkAgg")
    except Exception:
        pass

    df = pd.read_csv(csv_path)
    required = list(REQUIRED_BASE)
    if mode == "mar":
        required.append(MAR_Y_COL)
    missing = [c for c in required if c not in df.columns]
    if missing:
        print(f"ERROR: CSV missing required columns for {mode} mode: {missing}")
        if mode == "mar":
            print("Hint: run mass_accumulation_rate.py first to produce a MassDepth CSV.")
        return 3

    core_name = args.core or csv_path.stem
    _editor = CurveEditor(df, mode, core_name, csv_path.parent)
    plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
