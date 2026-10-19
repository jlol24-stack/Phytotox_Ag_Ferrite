# Phytotoxicity Analysis Pipeline for Ag@CoFe2O4 Nanoparticles in Wheat Plants

An advanced ecotoxicological data data-analytics pipeline implemented in Python (Jupyter Notebook) for processing phytotoxicity, bioaccumulation, and multivariate physiological biomarker data from crop assays exposed to core-shell nanoparticles.

## 📝 Citation Policy & Academic Integrity

If you utilize this notebook, code snippets, or analytical pipeline in your academic or professional research, **you are required to cite both this source code repository and the original research paper** upon which this mathematical implementation is built:

### 1.López‑Luna, J., Cruz-Fernández, S., Vázquez-Hipólito, V., Martinez‑Vargas S., Martínez, I.A., Mijangos‑Ricardez, O.F., Solís‑Domínguez, F.A., Cuevas‑Díaz, M.C. 2026. Phytotoxicity Analysis Pipeline for Ag@CoFe2O4 Nanoparticles in Wheat Plants. Version 1.0.0 10.5281/zenodo.23377918.

### 2. Original Research Article (Mandatory Reference)
> **Phytotoxicity and upper localization of Ag@CoFe2O4 nanoparticles in wheat plants**
> Jaime López-Luna, Soledad Cruz-Fernández, Donald Stewart Mills, Arturo Isaías Martínez-Enríquez, Fernando Amilcar Solís-Domínguez, María del Carmen Ángeles González-Chávez, Rogelio Carrillo-González, Sergio Martinez-Vargas, Oscar Francisco Mijangos-Ricardez, & María del Carmen Cuevas-Díaz.
> *Environmental Science and Pollution Research* (2020) 27:1923–1940.
> **DOI:** [10.1007/s11356-019-06668-9]

---

## 🔬 Computational Workflow

The notebook structures the analytical pipeline into 4 distinct experimental processing stages:

1. **Stage 1: Data Loading & Matrix LOD Treatment**
   * Imports raw matrix observation spreadsheets (`phytotox_ferrite.xlsx`).
   * Automatically handles data-censoring by substituting values below the ICP-OES instrument limit of detection using the robust scientific **\(\text{LOD}/2\)** imputation rule.
   * Segregates data arrays by vegetative organ compartments (`shoot` vs. `root`).

2. **Stage 2: Bioaccumulation Profiling**
   * Groups and compiles elemental translocations and total burdens (`Fe`, `Co`, `Ag`) across the explicit application concentration steps.

3. **Stage 3: Non-Linear Dose-Response Regression**
   * Fits a **3-Parameter Log-Logistic Model (LL.3)** to establish the median inhibitory concentration (\(IC_{50}\)) for elongation.
   * Integrates an automated **Dunnett's Post-Hoc Test** (\(p < 0.05\)) against the unamended control baseline to explicitly isolate the **NOAEC** and **LOAEC** regulatory safety thresholds.

4. **Stage 4: Multivariate Physiological System Dynamics**
   * Computes Pearson cross-correlation matrices connecting metal burdens to enzymatic stress networks.
   * Performs an **algebraic Singular Value Decomposition (SVD)** for dimensional reduction.
   * Generates custom publication-quality **PCA Biplots** embedded with automated axis expansion guards to guarantee boundary text labels (`Growth`, `CAT`, `Car`) remain safely contained inside the visible plot frame.

## ⚙️ Mathematical & Interpretation Notes

* **Hormesis Constraints:** The non-linear regression model (LL.3) focuses exclusively on monotonic inhibitory decay. Low-dose biostimulatory/hormetic trends are treated as variance rather than a separate biphasic sub-fit.
* **PCA Vector Orientations:** Axis orientations are subject to the standard "Sign Flipping" phenomenon depending on row formatting and library versions. Both mirror projections remain 100% mathematically correct; the underlying biological interactions, vector angles, and sample distances are perfectly invariant.

## 🛠️ Requirements & Environment Setup

Ensure the following core Python ecosystem modules are active in your Jupyter/Cursor kernel:

```bash
pip install pandas numpy matplotlib scipy openpyxl
```
