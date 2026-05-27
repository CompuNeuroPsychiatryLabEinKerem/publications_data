# Spatiotemporal Disorientation and Alzheimer's Disease Risk Factors

Code accompanying the manuscript:

> **Spatiotemporal disorientation in cognitively unimpaired adults is related to non-vascular risk factors of Alzheimer's disease**
>
> Gazit L, Yundler E, Elias U, Marshall GA & Arzy S

## Overview

Preprocessing and analysis code for the study. We tested whether modifiable risk factors for Alzheimer's disease (per the 2024 Lancet Commission) are associated with orientation performance in cognitively unimpaired adults assessed via the Clara mobile platform.

## Repository Structure

```
preprocessing/
    preprocessing.py          Raw data → merged participant dataset
    efficacy_scoring.py       Efficacy Score (ES) computation per game session

analysis/
    risk_factor_definitions.py    Binary risk-factor coding (13 Lancet factors + air pollution)
    propensity_score_analysis.py  Overlap-weighted doubly robust regression
    air_pollution_scoring.py      Comprehensive Pollution Index (CPI) from monitoring-station data

requirements.txt              Python dependencies
```

## Installation

```bash
pip install -r requirements.txt
```

## Running the Pipeline

Three stages, run in order. All scripts accept `--help`.

### Step 1: Preprocess raw data

Reads the raw Excel export and merges clients, games, and medical questionnaire data.

```bash
cd preprocessing/
python preprocessing.py --input /path/to/ShomeroshTables.xlsx --output-dir .
```

**Input:** `ShomeroshTables.xlsx` (Clients, Games, Medical sheets)
**Output:** `merged_data.pkl`

### Step 2: Compute Efficacy Scores

Clips response times to [0.3, 12] seconds, computes ES = success_rate / response_time per game session, and retains only the first session per participant per game type.

```bash
cd preprocessing/
python efficacy_scoring.py --input merged_data.pkl --output-dir .
```

**Input:** `merged_data.pkl` (from Step 1)
**Output:** `merged_df_ES.pkl`, `merged_df_ES.csv`

### Step 3a: Air pollution scoring (optional)

Scores Israeli monitoring stations using a Comprehensive Pollution Index (CPI) normalised to WHO 2021 guidelines. Only needed if air pollution is used as a risk factor.

```bash
cd analysis/
python air_pollution_scoring.py --data /path/to/air_quality_general_stations_2017.csv --output-dir .
```

**Input:** Monitoring-station CSV with annual pollutant concentrations
**Output:** `station_scores.csv`, `area_pollution_summary.csv`

### Step 3b: Propensity score analysis

Runs overlap-weighted doubly robust regression for all risk factors × cognitive measures.

```bash
cd analysis/
python propensity_score_analysis.py \
    --data ../preprocessing/merged_df_ES.csv \
    --output-dir results/
```

To include air pollution as a risk factor, add:
```bash
    --air-pollution-dir /path/to/dir/with/city_area_mapping.csv
```

**Input:** `merged_df_ES.csv` (from Step 2), optionally `city_area_mapping.csv`
**Output:** `results/propensity_results.csv` — one row per risk factor × measure, columns for both approaches

## Methods

### 1. Preprocessing

Raw data exported from the Clara app (Clients, Games, Medical sheets):

1. **Account filtering** — removed test accounts; kept ages 18–100.
2. **Medical questionnaire** — pivoted to one row per participant; merged gender-specific columns; computed BMI; clipped extreme sleep-hours values.
3. **Efficacy Scores** — ES = mean(success rate) / mean(response time) per game session. Response times clipped to 0.3–12 s. Kept only the first session per participant per game type.
4. **Exclusions** — dropped participants missing age or sex (needed for propensity model).

### 2. Risk Factor Definitions

Thirteen binary risk factors from self-reported questionnaire data, following the Lancet Commission (2024):

| Category | Risk Factors |
|---|---|
| **Non-vascular** | Low education (≤ 12 yrs), hearing loss, depression, anxiety, sleep disorders, physical inactivity, social isolation, excessive alcohol |
| **Vascular** | Hypertension, diabetes, high cholesterol, smoking (AHE), obesity (BMI ≥ 30) |

Air pollution was scored using a Comprehensive Pollution Index (CPI) from Israeli monitoring stations normalised to WHO 2021 guidelines, mapped to participants by city of residence. CPI > 2.0 = very high exposure.

### 3. Statistical Analysis

We used **overlap-weighted doubly robust regression**:

1. Propensity score: PS = P(Exposure | Age, Age², Sex, Age × Sex) via logistic regression.
2. Overlap weights: w = (1 − PS) for exposed, w = PS for unexposed.
3. Compared nested WLS models:
   - Null:  Score ~ Age + Age² + Sex + Age × Sex
   - Full:  Score ~ Exposure + Age + Age² + Sex + Age × Sex
4. Effect size = incremental variance explained (ΔR²).

Repeated for all risk factors × 7 cognitive measures (4 individual ES domains + 3 composites). Primary outcome: spatiotemporal orientation (mean of Places and Events).

## Data Availability

Data is not publicly available for privacy reasons but may be shared for research purposes upon request.

## References

- Livingston G, Huntley J, Liu KY, et al. Dementia prevention, intervention, and care: 2024 report of the Lancet standing Commission. *Lancet* 2024; 404: 572–628.
- Li F, Morgan KL, Zaslavsky AM. Balancing covariates via propensity score weighting. *JASA* 2018; 113(521): 390–400.
- WHO Global Air Quality Guidelines, 2021.
