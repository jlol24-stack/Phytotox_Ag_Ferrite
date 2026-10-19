# ============================================================
# PHYTOTOXICITY ANALYSIS OF WHEAT SEEDLINGS EXPOSED TO SILVER COATED COBALT FERRITE NPs
# File: phytotox_ferrite.xlsx
# Author: Jaime López-Luna
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import dunnett
import os

# Configure Matplotlib fonts and enable LaTeX rendering support globally
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rcParams['mathtext.fontset'] = 'dejavusans'  # Enables clear subscripts in labels

# =====================================================================
# STAGE 1: DATA LOADING, MATRIX LOD TREATMENT (µg/g) & SEGREGATION
# =====================================================================
file_path = 'phytotox_ferrite.xlsx'

if not os.path.exists(file_path):
    print(f"❌ Error: The file '{file_path}' was not found in: {os.getcwd()}")
    exit()

try:
    df = pd.read_excel(file_path)
    print("✅ Stage 1: Data successfully loaded into Cursor.")
except Exception as e:
    print(f"❌ Error loading file: {e}")
    exit()

metals = ['Fe', 'Co', 'Ag']
print("\n--- ICP-OES Matrix LOD Treatment (µg/g) ---")
for metal in metals:
    df[metal] = pd.to_numeric(df[metal], errors='coerce')
    min_quantified = df[df[metal] > 0][metal].min()
    operational_lod_half = min_quantified / 2
    under_lod_condition = df['Growth'].notna() & (df[metal].isna() | (df[metal] == 0))
    df.loc[under_lod_condition, metal] = operational_lod_half
    print(f"🔹 {metal:2s} -> Imputed (LOD/2): {operational_lod_half:8.4f} µg/g")

df_shoot = df[df['Organ'] == 'shoot'].copy()
df_root = df[df['Organ'] == 'root'].copy()

# =====================================================================
# STAGE 2: BIOACCUMULATION PROFILE ANALYSIS
# =====================================================================
def analyze_bioaccumulation(df_organ, organ_name):
    print(f"\n==================================================")
    print(f"STAGE 2: BIOACCUMULATION IN {organ_name.upper()} (µg/g)")
    print(f"==================================================")
    summary = df_organ.groupby('Dose')[metals].mean()
    print(summary.round(3))

analyze_bioaccumulation(df_shoot, "shoot")
analyze_bioaccumulation(df_root, "root")

# =====================================================================
# STAGE 3: DOSE-RESPONSE MODELING & IC50 CALCULATION (LL.3 MODEL)
# =====================================================================
def log_logistic_3p(dose, top, ic50, hill_slope):
    """Standard Log-Logistic 3-Parameter Model (LL.3) - Bottom fixed at 0"""
    return top / (1 + (dose / ic50) ** hill_slope)

def fit_dose_response(df_organ, organ_name):
    print(f"\n==================================================")
    print(f"STAGE 3: TOXICOLOGICAL PROFILE - LL.3 & DUNNETT ({organ_name.upper()})")
    print(f"==================================================")
    
    df_fit = df_organ.dropna(subset=['Growth'])
    X = df_fit['Dose'].values
    Y = df_fit['Growth'].values
    
    # 3.1 DUNNETT'S TEST FOR NOAEC & LOAEC DETERMINATION
    control_growth = df_fit[df_fit['Dose'] == 0]['Growth'].values
    exposed_doses = sorted([d for d in df_fit['Dose'].unique() if d > 0])
    
    noaec_val = 0.0
    loaec_val = None
    doses_data = [df_fit[df_fit['Dose'] == dose]['Growth'].values for dose in exposed_doses]
    
    try:
        res_dunnett = dunnett(*doses_data, control=control_growth, alternative='less')
        p_values = res_dunnett.pvalue
        
        for idx, dose in enumerate(exposed_doses):
            p_val = p_values[idx]
            if p_val < 0.05:
                loaec_val = float(dose)
                if idx > 0:
                    noaec_val = float(exposed_doses[idx-1])
                break
        else:
            noaec_val = float(exposed_doses[-1])
            
        print(f"📋 Dunnett's Test Statistical Thresholds:")
        print(f"   • NOAEC (No Observed Adverse Effect Concentration): {noaec_val} mg/kg")
        print(f"   • LOAEC (Lowest Observed Adverse Effect Concentration): {f'{loaec_val} mg/kg' if loaec_val else 'No significant effect'} ⚠️")
        print("   -------------------------------------------------")
    except Exception as dun_err:
        print(f"ℹ️ Dunnett's test skipped or needs replication check: {dun_err}")

    # 3.2 LOG-LOGISTIC NON-LINEAR REGRESSION (IC50)
    p0 = [100.0, np.median(X), 1.0]
    try:
        popt, pcov = curve_fit(log_logistic_3p, X, Y, p0=p0, bounds=(0, [120, np.inf, np.inf]))
        top, ic50, slope = popt
        
        # Extracted safely as individual floats to fix string formatting conflicts
        perr = np.sqrt(np.diag(pcov))
        top_err = float(perr[0])
        ic50_err = float(perr[1])
        slope_err = float(perr[2])
        
        residuals = Y - log_logistic_3p(X, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((Y - np.mean(Y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        print(f"✅ Successful Mathematical Fit for {organ_name.capitalize()}:")
        print(f"   • Control Growth (Top)  : {top:.2f} ± {top_err:.2f} %")
        print(f"   • Phytotoxicity (IC50)  : {ic50:.2f} ± {ic50_err:.2f} mg/kg 🎯")
        print(f"   • Hill Slope Coefficient: {slope:.2f} ± {slope_err:.2f}")
        print(f"   • R-squared (R²)        : {r_squared:.4f}")
        
        # --- GRAPHICS GENERATION WITH GEOMETRIC DUNNETT LINES ---
        plt.figure(figsize=(5.5, 4.5))
        dot_color = 'darkgreen' if organ_name=='shoot' else 'saddlebrown'
        
        # Replicates scatter and curve fit
        plt.scatter(X, Y, color=dot_color, label='Replicates', alpha=0.7, edgecolors='black', zorder=3)
        x_curve = np.linspace(0, max(X), 500)
        plt.plot(x_curve, log_logistic_3p(x_curve, *popt), 'r-', linewidth=2, label=f'Fit ($IC_{{50}}$ = {ic50:.1f} mg/kg)', zorder=4)
        
        # Plot NOAEC regulatory vertical line (Green)
        plt.axvline(x=noaec_val, color='green', linestyle='--', linewidth=1.2, label=f'NOAEC ({noaec_val:.0f} mg/kg)', zorder=2)
        
        # Plot LOAEC regulatory vertical line (Red) if triggered
        if loaec_val is not None:
            plt.axvline(x=loaec_val, color='darkred', linestyle='--', linewidth=1.2, label=f'LOAEC ({loaec_val:.0f} mg/kg)', zorder=2)
        
        plt.title(f'Dose-Response (LL.3) - {organ_name.capitalize()}', fontsize=11, fontweight='bold')
        plt.xlabel(r'Concentration of $CoFe_2O_4@Ag$ ($mg/kg$)')
        plt.ylabel('Growth (% of Control)')
        plt.grid(True, linestyle=':', alpha=0.5, zorder=1)
        plt.legend(fontsize=9, loc='upper right')
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"❌ Curve fitting failed for {organ_name}: {e}")

fit_dose_response(df_shoot, "shoot")
fit_dose_response(df_root, "root")

# =====================================================================
# STAGE 4: PHYSIOLOGICAL RESPONSE (PCA BIPLOT ONLY - FIXED AXIS LIMITS)
# =====================================================================
def analyze_physiological_response(df_organ, organ_name):
    print(f"\n==================================================")
    print(f"STAGE 4: PHYSIOLOGICAL ANALYSIS ({organ_name.upper()})")
    print(f"==================================================")
    
    pigments = ['Chla', 'Chlb', 'Chltot', 'Car']
    biomarkers = ['Prot', 'CAT', 'APX', 'GPX', 'SOD']
    
    # Dynamic selection of parameters per organ matrix
    valid_pigments = [p for p in pigments if p in df_organ.columns and df_organ[p].notna().any()]
    valid_biomarkers = [b for b in biomarkers if b in df_organ.columns and df_organ[b].notna().any()]
    all_targets = valid_pigments + valid_biomarkers
    
    if not all_targets:
        print(f"No active physiological markers found for {organ_name}.")
        return
        
    corr_matrix = df_organ[metals + all_targets].corr()
    metal_vs_targets = corr_matrix.loc[metals, all_targets]
    
    print("Pearson Correlation Coefficients (r) against Metal Burdens:")
    print(metal_vs_targets.round(3))
    
    features = metals + ['Growth'] + all_targets
    df_sub = df_organ[features].dropna()
    data_matrix = df_sub.values
    
    mean_vec = np.mean(data_matrix, axis=0)
    std_vec = np.std(data_matrix, axis=0)
    std_vec[std_vec == 0] = 1.0
    standardized_data = (data_matrix - mean_vec) / std_vec
    
    U, S, Vt = np.linalg.svd(standardized_data, full_matrices=False)
    variance_explained = (S ** 2) / np.sum(S ** 2)
    
    pc1_var = float(variance_explained[0] * 100)
    pc2_var = float(variance_explained[1] * 100)
    
    scores = standardized_data @ Vt.T
    loadings = Vt.T
    
    plt.figure(figsize=(5.5, 4.5))
    doses = df_organ.loc[df_sub.index, 'Dose'].values
    scatter = plt.scatter(scores[:, 0], scores[:, 1], c=doses, cmap='viridis', edgecolors='black', s=40, zorder=3)
    plt.colorbar(scatter, label="Exposure Dose (mg/kg)")
    
    scale_factor = np.max(np.abs(scores[:, 0:2])) / np.max(np.abs(loadings[:, 0:2])) * 0.7
    
    # Track the extreme coordinates of all objects to adjust boundaries
    max_x, min_x = np.max(scores[:, 0]), np.min(scores[:, 0])
    max_y, min_y = np.max(scores[:, 1]), np.min(scores[:, 1])
    
    # Process and map vector projection loading paths
    for i, feature_name in enumerate(features):
        x_val = loadings[i, 0] * scale_factor
        y_val = loadings[i, 1] * scale_factor
        
        # Update extreme bounding limits for layout tracking
        max_x, min_x = max(max_x, x_val * 1.2), min(min_x, x_val * 1.2)
        max_y, min_y = max(max_y, y_val * 1.2), min(min_y, y_val * 1.2)
        
        plt.arrow(0, 0, x_val, y_val, color='darkred', head_width=0.12, head_length=0.12, linewidth=1.1, zorder=4)
        plt.text(x_val * 1.15, y_val * 1.15, feature_name, color='maroon', ha='center', va='center', fontweight='bold', fontsize=9)
        
    plt.title(f'PCA Biplot - {organ_name.capitalize()} Bio-System', fontsize=11, fontweight='bold')
    plt.xlabel(f'Principal Component 1 ({pc1_var:.1f}%)')
    plt.ylabel(f'Principal Component 2 ({pc2_var:.1f}%)')
    
    # FIX: Dynamically expand the view limits by 25% safety padding to accommodate text blocks
    plt.xlim(min_x * 1.25, max_x * 1.25)
    plt.ylim(min_y * 1.25, max_y * 1.25)
    
    plt.axhline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    plt.grid(True, linestyle=':', alpha=0.5, zorder=1)
    plt.tight_layout()
    plt.show()

# Parallel execution for both vegetative fractions
analyze_physiological_response(df_shoot, "shoot")
analyze_physiological_response(df_root, "root")

print("\n🚀 All 4 mathematical stages successfully processed with Dunnett lines.")








