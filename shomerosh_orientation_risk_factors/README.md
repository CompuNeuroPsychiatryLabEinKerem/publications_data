# Spatiotemporal Disorientation and Alzheimer's Disease Risk Factors

Code accompanying the manuscript:

> **Spatiotemporal disorientation in cognitively unimpaired adults is related to non-vascular risk factors of Alzheimer's disease**
>
> Gazit L, Yundler E, Elias U, Marshall GA & Arzy S

## Overview

This repository contains the data-processing pipeline and statistical analysis code used in the study. We examined associations between 14 modifiable risk factors identified by the Lancet Commission (2024) and orientation performance measured via the Clara/Shomerosh mobile platform in a large sample of cognitively unimpaired community-dwelling adults.

## Repository Structure

```
preprocessing/
    preprocessing.py          Raw data → merged participant dataset
    efficacy_scoring.py       Efficacy Score (ES) computation per game session

analysis/
    risk_factor_definitions.py    Binary risk-factor coding (13 Lancet factors + air pollution)
    propensity_score_analysis.py  Overlap-weighted doubly robust regression
    air_pollution_scoring.py      Comprehensive Pollution Index (CPI) from monitoring-station data
```

## Methods

### 1. Preprocessing

Participants were recruited through the Clara mobile app (locally adapted as *Shomerosh*; www.shomerosh.com). Raw data from the app export (Clients, Games, Medical sheets) was processed as follows:

1. **Account filtering** — test/developer accounts removed; ages restricted to 18–100.
2. **Medical questionnaire** — responses pivoted to one row per participant; gender-specific columns merged; BMI computed from validated height and weight; extreme sleep-hours values clipped.
3. **Efficacy Scores** — for each game session, ES = mean(success rate) / mean(response time). Response times were clipped to 0.3–12 seconds. Only the first session per participant per game type (Places, Events, People, Lexical) was retained.
4. **Exclusions** — participants missing age or sex were excluded, as these are required for the propensity score model (final analysis sample: N ≈ 2,908).

### 2. Risk Factor Definitions

Thirteen binary risk factors were created from self-reported medical and lifestyle questionnaire data following Lancet Commission (2024) definitions:

| Category | Risk Factors |
|---|---|
| **Non-vascular** | Low education (≤ 12 yrs), hearing loss, depression, anxiety, sleep disorders, physical inactivity, social isolation, excessive alcohol |
| **Vascular** | Hypertension, diabetes, high cholesterol, smoking (AHE), obesity (BMI ≥ 30) |

Air pollution exposure was derived from a Comprehensive Pollution Index (CPI) based on Israeli monitoring-station data, normalised against WHO Air Quality Guidelines (2021). Participants were assigned the CPI of their self-reported city of residence; CPI > 2.0 was classified as very high exposure.

### 3. Statistical Analysis

Associations between each risk factor and orientation performance were estimated using **overlap-weighted doubly robust regression**:

1. A **propensity score** PS = P(Exposure | Age, Age², Sex, Age × Sex) was estimated via logistic regression.
2. **Overlap weights** were computed: w = (1 − PS) for exposed, w = PS for unexposed participants.
3. Two nested **weighted least-squares** models were compared:
   - Null:  Score ~ Age + Age² + Sex + Age × Sex
   - Full:  Score ~ Exposure + Age + Age² + Sex + Age × Sex
4. The incremental variance explained (ΔR²) was used as the effect-size measure.

This procedure was repeated for all 14 risk factors × 7 cognitive measures (4 individual domains + 3 composite scores). Spatiotemporal orientation (mean of Places and Events) was the primary outcome.

## Data Availability

Data is not made publicly available for participants' privacy reasons and may be available for scientific research upon request, following verification.

## Requirements

```
pip install -r requirements.txt
```

## References

- Livingston G, Huntley J, Liu KY, et al. Dementia prevention, intervention, and care: 2024 report of the Lancet standing Commission. *Lancet* 2024; 404: 572–628.
- Li F, Morgan KL, Zaslavsky AM. Balancing covariates via propensity score weighting. *JASA* 2018; 113(521): 390–400.
- WHO Global Air Quality Guidelines, 2021.
