"""
Comprehensive Pollution Index (CPI) scoring.

Normalises pollutant concentrations from Israeli monitoring stations
against WHO Air Quality Guidelines (2021).  Pollutant weights follow
the 2024 Lancet Commission.

Pollutant weights:
    PM2.5  2.0   (primary dementia risk factor)
    PM10   1.5
    NO2    1.2
    O3     1.0
    SO2    0.8
    CO     0.5

CPI = sum(pollutant_index × weight) / sum(weights)

Risk classification:
    CPI <= 1.0   Low (meets WHO guidelines)
    1.0–2.0      Moderate
    2.0–3.0      Elevated
    3.0–5.0      High
    > 5.0        Very high

References:
    Livingston G et al. Lancet 2024; 404: 572-628.
    WHO Global Air Quality Guidelines, 2021.
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

WHO_GUIDELINES = {
    'PM2.5': {'annual_avg': 5.0,  'weight': 2.0},
    'PM10':  {'annual_avg': 15.0, 'weight': 1.5},
    'NO2':   {'annual_avg': 10.0, 'weight': 1.2},
    'O3':    {'annual_avg': 60.0, 'weight': 1.0},
    'SO2':   {'annual_avg': 20.0, 'weight': 0.8},
    'CO':    {'annual_avg': 1.0,  'weight': 0.5},   # mg/m³
}

COL_MAP = {
    'PM2.5': 'PM2.5_annual_avg_ug',
    'PM10':  'PM10_annual_avg_ug',
    'NO2':   'NO2_annual_avg_ug',
    'O3':    'O3_annual_avg_ug',
    'SO2':   'SO2_annual_avg_ug',
    'CO':    'CO_annual_avg_mg',
}

RISK_LEVELS = [
    ('Low',       0, 1.0),
    ('Moderate',  1.0, 2.0),
    ('Elevated',  2.0, 3.0),
    ('High',      3.0, 5.0),
    ('Very High', 5.0, float('inf')),
]


def compute_cpi(row):
    """Return CPI for one monitoring station."""
    w_sum, w_total = 0.0, 0.0
    for poll, guide in WHO_GUIDELINES.items():
        col = COL_MAP.get(poll)
        if col and col in row.index and pd.notna(row[col]):
            idx = row[col] / guide['annual_avg']
            w_sum += idx * guide['weight']
            w_total += guide['weight']
    return w_sum / w_total if w_total > 0 else np.nan


def classify(cpi):
    if pd.isna(cpi):
        return 'Unknown'
    for label, lo, hi in RISK_LEVELS:
        if lo < cpi <= hi:
            return label
    return 'Unknown'


def score_stations(csv_path):
    """Score all monitoring stations and return a summary DataFrame."""
    df = pd.read_csv(csv_path)
    df['CPI'] = df.apply(compute_cpi, axis=1)
    df['risk_level'] = df['CPI'].apply(classify)
    return df


def area_summary(scored_df):
    """Aggregate CPI to geographic-area level."""
    return scored_df.groupby('Area').agg(
        CPI_mean=('CPI', 'mean'),
        CPI_median=('CPI', 'median'),
        CPI_std=('CPI', 'std'),
        n_stations=('CPI', 'count'),
    ).sort_values('CPI_mean', ascending=False).round(2)


def main(csv_path, output_dir='.'):
    import os
    os.makedirs(output_dir, exist_ok=True)
    scored = score_stations(csv_path)
    scored.to_csv(os.path.join(output_dir, 'station_scores.csv'), index=False)
    summary = area_summary(scored)
    summary.to_csv(os.path.join(output_dir, 'area_pollution_summary.csv'))
    print(f"Scored {len(scored)} stations, {summary.shape[0]} areas")
    return scored, summary


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Score air-pollution monitoring stations')
    parser.add_argument('--data', default='air_quality_general_stations_2017.csv',
                        help='Path to monitoring-station CSV')
    parser.add_argument('--output-dir', default='.',
                        help='Directory for output files')
    args = parser.parse_args()
    main(args.data, args.output_dir)
