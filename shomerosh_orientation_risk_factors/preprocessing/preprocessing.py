"""
Data preprocessing for Shomerosh/Clara mobile platform.

Loads raw data (Clients, Games, Medical sheets), filters test accounts and
invalid ages (< 18 or > 100), pivots medical questionnaire responses,
computes BMI, and produces a merged per-user dataset (merged_data.pkl).
"""

import pandas as pd
import numpy as np
import os
import pickle
from datetime import date


def load_data(input_file):
    """Load Clients, Games, and Medical sheets from the raw Excel export."""
    clients = pd.read_excel(input_file, sheet_name='Clients')
    games = pd.read_excel(input_file, sheet_name='Games')
    medical = pd.read_excel(input_file, sheet_name='Medical')
    return clients, games, medical


def filter_clients(clients, date_of_interest=None):
    """Remove test accounts and participants with implausible ages."""
    df = clients[clients['testAccount'] != True].copy()
    df['updatedAt'] = pd.to_datetime(df['updatedAt'], format='%d.%m.%Y %H:%M:%S')

    if date_of_interest:
        df = df[df['updatedAt'] >= date_of_interest]

    df = df[(df['Age'].isna()) | ((df['Age'] >= 18) & (df['Age'] <= 100))]
    return df


def medical_preprocessing(medical):
    """Pivot medical questionnaire, synchronise gender-specific columns, compute BMI."""
    medical['Answer'] = medical['Answer'].apply(lambda x: f'{x:.0f}')
    medical = medical.dropna(axis=0)
    medical = medical.drop_duplicates(subset=['userID', 'Coding'], keep='first')
    medical = pd.pivot_table(medical, index='userID', columns='Coding',
                             values='Answer', aggfunc='first')

    # Combine gender-specific columns (e.g. smoke_perDay_f / smoke_perDay_m)
    factor_names = set(col[:-2] for col in medical.columns
                       if col.endswith('_f') or col.endswith('_m'))
    for factor in factor_names:
        f_col, m_col = f"{factor}_f", f"{factor}_m"
        if f_col in medical.columns and m_col in medical.columns:
            medical[factor] = medical[f_col].fillna(medical[m_col])
            medical.drop([f_col, m_col], axis=1, inplace=True)

    # BMI from validated height / weight
    medical['metrics_height_cm'] = pd.to_numeric(medical['metrics_height_cm'], errors='coerce')
    medical['metrics_weight_kg'] = pd.to_numeric(medical['metrics_weight_kg'], errors='coerce')
    medical['BMI'] = np.nan
    valid = medical['metrics_height_cm'].between(130, 220) & \
            medical['metrics_weight_kg'].between(30, 200)
    medical.loc[valid, 'BMI'] = (
        medical.loc[valid, 'metrics_weight_kg'] /
        (medical.loc[valid, 'metrics_height_cm'] / 100) ** 2
    ).round(2)

    # Clip extreme self-reported sleep hours
    medical['sleep_hours'] = pd.to_numeric(medical['sleep_hours'], errors='coerce')
    medical.loc[(medical['sleep_hours'] > 18) |
                (medical['sleep_hours'] < 2), 'sleep_hours'] = np.nan

    return medical


def merge_all_dfs(df_clients, df_games, df_medical):
    """Inner-join games with filtered clients; left-join medical data."""
    merged = pd.merge(df_games, df_clients, on='userID', how='inner')
    merged = pd.merge(df_medical, merged, on='userID', how='right')
    return merged


def main(input_file='ShomeroshTables.xlsx', output_dir='.'):
    os.makedirs(output_dir, exist_ok=True)
    clients, games, medical = load_data(input_file)

    clients = filter_clients(clients)
    medical = medical_preprocessing(medical)
    merged_df = merge_all_dfs(clients, games, medical)

    out_path = os.path.join(output_dir, 'merged_data.pkl')
    with open(out_path, 'wb') as f:
        pickle.dump(merged_df, f)

    print(f"Saved merged_data.pkl — {merged_df['userID'].nunique()} unique users")
    return merged_df


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Preprocess raw Shomerosh data')
    parser.add_argument('--input', default='ShomeroshTables.xlsx',
                        help='Path to raw Excel export')
    parser.add_argument('--output-dir', default='.',
                        help='Directory for output files')
    args = parser.parse_args()
    main(args.input, args.output_dir)
