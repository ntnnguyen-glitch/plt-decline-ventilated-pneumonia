# Early platelet decline and 28-day mortality in mechanically ventilated patients with pneumonia

R code for a retrospective cohort study (landmark analysis at 96 h after ICU admission) using MIMIC-IV v3.1 (development) and eICU-CRD v2.0 (external validation). No machine-learning models are used. Reporting follows STROBE.

## Data access
Both databases require PhysioNet credentialed access and a signed data use agreement. No patient-level data are included in this repository.

## Requirements
R ≥ 4.3 (analyses were run in R 4.6.1) with: data.table, DBI, duckdb, survival, splines, mice, quantreg, WeightIt, MatchIt, pROC, EValue, sandwich, lmtest, ggplot2, patchwork, ragg, officer, flextable.

## Run
```bash
MIMIC_DIR=/path/to/mimic-iv-3.1 EICU_DIR=/path/to/eicu-crd-2.0 ./run_all.sh
```

| Script | Content |
|---|---|
| `R/00_extract.R` | Extracts the required MIMIC-IV and eICU-CRD tables from raw csv.gz files into a local DuckDB file |
| `R/00b_extract_hemostasis.R` | Platelet count, INR and aPTT at all time points, anticoagulant and antiplatelet orders (with dose), blood-product transfusions, vasopressors, renal replacement therapy |
| `R/01_cohort.R` | Adults, first ICU stay, pneumonia, invasive ventilation started −6 to 24 h, alive at 96 h, baseline platelet (first value −6 to 24 h) and ≥1 value at 48–96 h; excludes haematological malignancy, ECMO and therapeutic anticoagulation within ±24 h (including LMWH/fondaparinux above prophylactic doses). Exposure: relative platelet decline (baseline − nadir 48–96 h)/baseline ≥30% |
| `R/02_analysis.R` | Multiple imputation by chained equations (20 datasets, Rubin's rules); 28-day mortality (logistic Models 1–3, Cox from the landmark), continuous/categorical decline, Kaplan–Meier, restricted cubic splines, incremental value; secondary outcomes (in-hospital, 90-day and 1-year mortality, ventilator-free days and ICU stay by median regression, platelet and red-cell transfusion); IPTW, doubly robust, overlap weights, propensity score matching, E-value, hospital-clustered SE; sensitivity analyses (thresholds, windows, time-varying organ support 24–96 h, pneumonia definition, exclusions, median imputation, complete-case) and subgroup analyses with interaction tests |
| `R/03_figures.R` | Figures 1–4 (TIFF, 300 dpi) |
| `R/04_tables.R` | Tables 1–3 and supplementary tables |

All outputs are written to `output/`.
