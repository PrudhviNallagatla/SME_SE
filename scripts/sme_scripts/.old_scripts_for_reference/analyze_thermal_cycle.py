#!/usr/bin/env python3
"""
analyze_thermal_cycle_v2.py

Analyzes the 'thermo_detailed.txt' output files.

VERSION 2: Includes automated sanity checks to validate the
           transformation temperatures without visual inspection.

This script:
1.  Finds all 'thermo_detailed.txt' files.
2.  Parses the data, separating cooling and heating ramps.
3.  Uses the derivative method (d(PE)/dT) to find (Ms, Mf, As, Af).
4.  ### NEW ### Runs programmatic sanity checks on the results.
5.  Saves a 4-panel plot for auditing/debugging.
6.  Saves a 'transformation_summary.txt' with the key values AND
    a clear success/failure flag.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import savgol_filter
import sys
import os

# ============================================
# CONFIGURATION
# ============================================
PROJECT_ROOT = Path("/home/rimuru/workspace")
INPUT_DATA_ROOT = PROJECT_ROOT / "data" / "raw_output"
ANALYSIS_OUTPUT_ROOT = PROJECT_ROOT / "data" / "analysis_output"

ANALYSIS_CONFIG = {
    # Savitzky-Golay filter window length (must be odd)
    # Larger window = more smoothing (good for noisy data)
    "WINDOW_LENGTH": 51,
    # Savitzky-Golay filter polynomial order (usually 2 or 3)
    "POLYORDER": 3,
    # Percentage of peak height to define start/end of transformation
    # Smaller value = wider transformation range
    "THRESHOLD_PERCENT": 0.05,  # 5%
    
    ### NEW ###
    # Minimum signal-to-noise ratio for the derivative peak.
    # A value > 5 is usually good. If you get warnings,
    # try increasing WINDOW_LENGTH or lowering this threshold.
    "MIN_PEAK_SIGNIFICANCE": 5.0 
}
# ============================================

# (The load_data and split_cycles functions are unchanged)
def load_data(filepath: Path) -> pd.DataFrame:
    """Loads the thermo_detailed.txt file into a pandas DataFrame."""
    try:
        header_line = ""
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    header_line = line.strip().lstrip("# ").split()
                    break
        if not header_line:
            raise ValueError("Could not parse header.")
            
        df = pd.read_csv(
            filepath,
            delim_whitespace=True,
            comment="#",
            names=header_line,
            skiprows=1
        )
        df = df.apply(pd.to_numeric, errors='coerce')
        return df.dropna()
        
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error loading data from {filepath}: {e}", file=sys.stderr)
        return None

def split_cycles(df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    """
    Splits the full dataframe into cooling and heating ramps.
    Filters out equilibration steps by looking at the temperature derivative.
    """
    df['dT'] = df['temp'].diff()
    # A more robust filter to find the main ramp
    cooling_df = df[df['dT'] < -1e-6].copy()
    heating_df = df[df['dT'] > 1e-6].copy()
    return cooling_df, heating_df


def analyze_transformation(
    cycle_df: pd.DataFrame, 
    property_col: str = 'pe', 
    temp_col: str = 'temp'
) -> (float, float, float, dict):
    """
    Analyzes a single ramp (cooling or heating) to find
    start, peak, and finish temperatures (e.g., Ms, M_peak, Mf).
    
    Returns: (T_start, T_peak, T_finish, analysis_data_dict)
    """
    cfg = ANALYSIS_CONFIG
    
    cycle_df = cycle_df.sort_values(by=temp_col)
    T = cycle_df[temp_col].values
    P = cycle_df[property_col].values
    
    if len(T) < cfg["WINDOW_LENGTH"]:
        print(
            f"  [WARN] Not enough data points ({len(T)}) for "
            f"smoothing window ({cfg['WINDOW_LENGTH']}). Skipping.",
            file=sys.stderr
        )
        return np.nan, np.nan, np.nan, {"failed": True}

    try:
        P_smooth = savgol_filter(
            P, 
            window_length=cfg["WINDOW_LENGTH"], 
            polyorder=cfg["POLYORDER"]
        )
    except ValueError as e:
        print(f"  [WARN] Error during smoothing: {e}. Skipping.", file=sys.stderr)
        return np.nan, np.nan, np.nan, {"failed": True}

    dPdT = np.gradient(P_smooth, T)
    
    dPdT_smooth = savgol_filter(
        dPdT, 
        window_length=cfg["WINDOW_LENGTH"], 
        polyorder=cfg["POLYORDER"]
    )
    
    # 5. Find the transformation peak (positive peak for both)
    peak_idx = np.argmax(dPdT_smooth)
    T_peak = T[peak_idx]
    
    # 6. Find start and finish temperatures
    baseline_val = np.min(dPdT_smooth)
    peak_height = dPdT_smooth[peak_idx] - baseline_val
    threshold = baseline_val + cfg["THRESHOLD_PERCENT"] * peak_height
    
    T_start, T_finish = np.nan, np.nan
    try:
        active_indices = np.where(dPdT_smooth > threshold)[0]
        if len(active_indices) == 0:
            raise ValueError("No transformation peak found above threshold.")
            
        T_start = T[active_indices[0]]
        T_finish = T[active_indices[-1]]
        
    except Exception as e:
        print(f"  [WARN] Could not find T_start/T_finish: {e}", file=sys.stderr)

    ### NEW: Calculate metrics for sanity checks ###
    # We define "noise" as the standard deviation of the *un-smoothed*
    # derivative, which is a good measure of the data's "jiggle".
    # A more robust noise measure: std of the residual
    noise_level = np.std(dPdT - dPdT_smooth) 
    peak_significance = peak_height / noise_level if noise_level > 1e-9 else 0.0

    analysis_data = {
        "T": T,
        "P_smooth": P_smooth,
        "dPdT_smooth": dPdT_smooth,
        "threshold": threshold,
        "peak_significance": peak_significance,
        "failed": False
    }
    
    return T_start, T_peak, T_finish, analysis_data


### NEW: Sanity Check Functions ###
def run_sanity_checks(analysis: dict, errors: list) -> bool:
    """
    Runs a series of automated checks on the analysis results.
    Returns True if OK, False if any check fails.
    Modifies the 'errors' list in-place.
    """
    cfg = ANALYSIS_CONFIG
    Ms, M_peak, Mf = analysis['cooling']['temps']
    As, A_peak, Af = analysis['heating']['temps']
    cool_data = analysis['cooling']['plot_data']
    heat_data = analysis['heating']['plot_data']
    
    is_ok = True

    # Check 1: Did the analysis run at all?
    if cool_data.get("failed", False) or heat_data.get("failed", False):
        errors.append("Analysis routine failed (e.g., not enough data).")
        is_ok = False
        # Stop here, other checks will fail
        return is_ok

    # Check 2: Were all temperatures found?
    temps = {'Ms': Ms, 'Mf': Mf, 'As': As, 'Af': Af}
    for name, T in temps.items():
        if np.isnan(T):
            errors.append(f"Could not determine '{name}' (result is NaN).")
            is_ok = False
            
    if not is_ok:
        # Stop here, physical ordering checks will fail
        return is_ok

    # Check 3: Physical Ordering
    if Ms <= Mf:
        errors.append(f"Physicality Error: Ms ({Ms:.1f} K) must be > Mf ({Mf:.1f} K).")
        is_ok = False
    if Af <= As:
        errors.append(f"Physicality Error: Af ({Af:.1f} K) must be > As ({As:.1f} K).")
        is_ok = False
    if As <= Mf:
        errors.append(f"Physicality Error: As ({As:.1f} K) must be > Mf ({Mf:.1f} K).")
        is_ok = False
        
    # Check 4: Peak Significance (Signal-to-Noise)
    cool_sig = cool_data.get('peak_significance', 0)
    heat_sig = heat_data.get('peak_significance', 0)
    
    if cool_sig < cfg["MIN_PEAK_SIGNIFICANCE"]:
        errors.append(
            f"Cooling peak is not significant (S/N: {cool_sig:.1f} "
            f"< threshold: {cfg['MIN_PEAK_SIGNIFICANCE']}). "
            "Data may be too noisy or smoothing is wrong."
        )
        is_ok = False
        
    if heat_sig < cfg["MIN_PEAK_SIGNIFICANCE"]:
        errors.append(
            f"Heating peak is not significant (S/N: {heat_sig:.1f} "
            f"< threshold: {cfg['MIN_PEAK_SIGNIFICANCE']}). "
            "Data may be too noisy or smoothing is wrong."
        )
        is_ok = False
        
    return is_ok

# (The plot_results function is unchanged)
def plot_results(
    cooling_df: pd.DataFrame, 
    heating_df: pd.DataFrame,
    analysis: dict,
    output_path: Path
):
    """Generates and saves a 4-panel analysis plot."""
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(f"Thermal Transformation Analysis: {output_path.parts[-3]}", fontsize=18)
    
    Ms, M_peak, Mf = analysis['cooling']['temps']
    As, A_peak, Af = analysis['heating']['temps']
    cool_data = analysis['cooling']['plot_data']
    heat_data = analysis['heating']['plot_data']

    # --- Plot 1: PE Hysteresis ---
    ax = axes[0, 0]
    ax.plot(cooling_df['temp'], cooling_df['pe'], '.-', label='Cooling (A->M)', color='blue', alpha=0.7)
    ax.plot(heating_df['temp'], heating_df['pe'], '.-', label='Heating (M->A)', color='red', alpha=0.7)
    if not np.isnan(Ms):
        ax.axvline(Ms, ls='--', color='blue', label=f'Ms: {Ms:.1f} K')
        ax.axvline(Mf, ls=':', color='blue', label=f'Mf: {Mf:.1f} K')
    if not np.isnan(As):
        ax.axvline(As, ls='--', color='red', label=f'As: {As:.1f} K')
        ax.axvline(Af, ls=':', color='red', label=f'Af: {Af:.1f} K')
    ax.set_xlabel("Temperature (K)", fontsize=12)
    ax.set_ylabel("Potential Energy (eV)", fontsize=12)
    ax.set_title("Potential Energy Hysteresis", fontsize=14)
    ax.legend()

    # --- Plot 2: Volume Hysteresis ---
    ax = axes[0, 1]
    ax.plot(cooling_df['temp'], cooling_df['vol'], '.-', label='Cooling (A->M)', color='blue', alpha=0.7)
    ax.plot(heating_df['temp'], heating_df['vol'], '.-', label='Heating (M->A)', color='red', alpha=0.7)
    ax.set_xlabel("Temperature (K)", fontsize=12)
    ax.set_ylabel("Volume (Å³)", fontsize=12)
    ax.set_title("Volume Hysteresis", fontsize=14)
    ax.legend()

    # --- Plot 3: Cooling Derivative (d(PE)/dT) ---
    ax = axes[1, 0]
    if 'T' in cool_data:
        ax.plot(cool_data['T'], cool_data['dPdT_smooth'], '-', color='blue', label='d(PE)/dT (smoothed)')
        ax.axhline(cool_data['threshold'], ls=':', color='gray', label=f'{ANALYSIS_CONFIG["THRESHOLD_PERCENT"]*100:.0f}% Threshold')
        if not np.isnan(Ms):
            ax.axvline(Ms, ls='--', color='k', label=f'Ms: {Ms:.1f} K')
            ax.axvline(M_peak, ls='-', color='gray', label=f'Peak: {M_peak:.1f} K')
            ax.axvline(Mf, ls=':', color='k', label=f'Mf: {Mf:.1f} K')
        # Add S/N to plot
        sig = cool_data.get('peak_significance', 0)
        ax.text(0.05, 0.95, f"Peak S/N: {sig:.1f}", transform=ax.transAxes, 
                ha='left', va='top', bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.5))
    ax.set_xlabel("Temperature (K)", fontsize=12)
    ax.set_ylabel("d(PE)/dT (eV/K)", fontsize=12)
    ax.set_title("Cooling Analysis (A → M)", fontsize=14)
    ax.legend()

    # --- Plot 4: Heating Derivative (d(PE)/dT) ---
    ax = axes[1, 1]
    if 'T' in heat_data:
        ax.plot(heat_data['T'], heat_data['dPdT_smooth'], '-', color='red', label='d(PE)/dT (smoothed)')
        ax.axhline(heat_data['threshold'], ls=':', color='gray', label=f'{ANALYSIS_CONFIG["THRESHOLD_PERCENT"]*100:.0f}% Threshold')
        if not np.isnan(As):
            ax.axvline(As, ls='--', color='k', label=f'As: {As:.1f} K')
            ax.axvline(A_peak, ls='-', color='gray', label=f'Peak: {A_peak:.1f} K')
            ax.axvline(Af, ls=':', color='k', label=f'Af: {Af:.1f} K')
        # Add S/N to plot
        sig = heat_data.get('peak_significance', 0)
        ax.text(0.05, 0.95, f"Peak S/N: {sig:.1f}", transform=ax.transAxes, 
                ha='left', va='top', bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.5))
    ax.set_xlabel("Temperature (K)", fontsize=12)
    ax.set_ylabel("d(PE)/dT (eV/K)", fontsize=12)
    ax.set_title("Heating Analysis (M → A)", fontsize=14)
    ax.legend()

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    save_path = output_path / "transformation_analysis.png"
    plt.savefig(save_path, dpi=200)
    plt.close(fig)
    print(f"    ✓ Saved analysis plot to {save_path}")


def save_summary(analysis: dict, analysis_ok: bool, errors: list, output_path: Path):
    """Saves the key temperature values and analysis status to a text file."""
    
    Ms, M_peak, Mf = analysis['cooling']['temps']
    As, A_peak, Af = analysis['heating']['temps']
    
    hysteresis_width = (Af - Ms) if not (np.isnan(Af) or np.isnan(Ms)) else np.nan
    hysteresis_center_shift = (A_peak - M_peak) if not (np.isnan(A_peak) or np.isnan(M_peak)) else np.nan

    summary_file = output_path / "transformation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"# Thermal Transformation Analysis Summary\n")
        f.write(f"# Structure: {output_path.parts[-3]}\n")
        f.write("-" * 60 + "\n")
        
        ### NEW: Analysis Status ###
        if analysis_ok:
            f.write("ANALYSIS STATUS: [✓ OK] Automated checks passed.\n")
        else:
            f.write("ANALYSIS STATUS: [✗ FAILED] Automated checks failed.\n")
            f.write("    REASONS:\n")
            for i, err in enumerate(errors, 1):
                f.write(f"      {i}. {err}\n")
            f.write("\n    WARNING: Do NOT trust these values. "
                    "Check 'transformation_analysis.png' and "
                    "tune ANALYSIS_CONFIG parameters in the script.\n")
                    
        f.write("-" * 60 + "\n\n")
        
        f.write(f"MARTENSITE TRANSFORMATION (COOLING)\n")
        f.write(f"  Ms (Start):   {Ms:.2f} K\n")
        f.write(f"  M_peak (Peak):  {M_peak:.2f} K\n")
        f.write(f"  Mf (Finish):  {Mf:.2f} K\n")
        f.write("\n")
        f.write(f"AUSTENITE TRANSFORMATION (HEATING)\n")
        f.write(f"  As (Start):   {As:.2f} K\n")
        f.write(f"  A_peak (Peak):  {A_peak:.2f} K\n")
        f.write(f"  Af (Finish):  {Af:.2f} K\n")
        f.write("\n")
        f.write(f"HYSTERESIS PROPERTIES\n")
        f.write(f"  Width (Af - Ms):         {hysteresis_width:.2f} K\n")
        f.write(f"  Peak Shift (A_peak - M_peak): {hysteresis_center_shift:.2f} K\n")
        f.write("-" * 60 + "\n")
        
        if analysis_ok:
            f.write("\nACTION REQUIRED:\n")
            f.write("  Update 'sme_DEFORMATION_FIXED.lmp' with these values:\n")
            # Provide a 20K safety margin
            f.write(f"  variable T_low  equal {Mf - 20.0:.1f}  # (Value below Mf: {Mf:.1f} K)\n")
            f.write(f"  variable T_high equal {Af + 20.0:.1f}  # (Value above Af: {Af:.1f} K)\n")
        else:
            f.write("\nACTION REQUIRED:\n")
            f.write("  Analysis failed. You must manually inspect the .png plot\n")
            f.write("  and tune the 'ANALYSIS_CONFIG' parameters at the top\n")
            f.write("  of this python script (analyze_thermal_cycle_v2.py).\n")
            
    print(f"    ✓ Saved summary to {summary_file}")


def main():
    """ Main driver function. """
    print("==========================================")
    print("Running Thermal Cycle Analysis (v2)")
    print(f"Searching for files in: {INPUT_DATA_ROOT}")
    print("==========================================")
    
    ANALYSIS_OUTPUT_ROOT.mkdir(exist_ok=True)
    thermo_files = list(INPUT_DATA_ROOT.glob("*/thermal_cycle/thermo_detailed.txt"))
    
    if not thermo_files:
        print(f"Error: No 'thermo_detailed.txt' files found.", file=sys.stderr)
        return
        
    print(f"Found {len(thermo_files)} files to analyze...")
    
    overall_success = True
    for thermo_file in thermo_files:
        structure_name = thermo_file.parts[-3]
        print(f"\nProcessing: {structure_name}")
        
        output_dir = ANALYSIS_OUTPUT_ROOT / structure_name / "thermal_cycle"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        df = load_data(thermo_file)
        if df is None:
            overall_success = False
            continue
            
        cooling_df, heating_df = split_cycles(df)
        if cooling_df.empty or heating_df.empty:
            print("  [ERROR] Could not find both cooling and heating ramps. Skipping.", file=sys.stderr)
            overall_success = False
            continue
            
        print("  Analyzing cooling ramp (A->M)...")
        Ms, M_peak, Mf, cool_data = analyze_transformation(cooling_df)
        
        print("  Analyzing heating ramp (M->A)...")
        As, A_peak, Af, heat_data = analyze_transformation(heating_df)
        
        analysis_results = {
            "cooling": {"temps": (Ms, M_peak, Mf), "plot_data": cool_data},
            "heating": {"temps": (As, A_peak, Af), "plot_data": heat_data}
        }
        
        # --- 5. ### NEW ### Run Sanity Checks ---
        print("  Running automated sanity checks...")
        analysis_errors = []
        analysis_ok = run_sanity_checks(analysis_results, analysis_errors)
        
        if analysis_ok:
            print("    [✓ ANALYSIS OK] All checks passed.")
        else:
            print(f"    [✗ ANALYSIS FAILED] {len(analysis_errors)} error(s) found:")
            for err in analysis_errors:
                print(f"      - {err}")
            overall_success = False

        # --- 6. Plot Results (always, for debugging) ---
        plot_results(cooling_df, heating_df, analysis_results, output_dir)
        
        # --- 7. Save Summary ---
        save_summary(analysis_results, analysis_ok, analysis_errors, output_dir)

    print("\n==========================================")
    if overall_success:
        print("✓ Analysis complete. All structures passed checks.")
    else:
        print("! Analysis complete. One or more structures FAILED automated checks.")
        print("  Please review the 'transformation_summary.txt' and '.png' files")
        print(f"  in {ANALYSIS_OUTPUT_ROOT} for the failed structures.")
    print("==========================================")


if __name__ == "__main__":
    main()