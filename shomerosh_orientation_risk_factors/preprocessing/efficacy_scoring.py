"""
Efficacy Score (ES) calculation per game session.

Requires: merged_data.pkl (output of preprocessing.py).

For each game session the Efficacy Score is defined as:
    ES = mean(success_rate) / mean(response_time)

Response times are clipped to [0.3, 12] seconds before scoring.
Only the first session per participant per game type is retained.
Output: merged_df_ES.pkl (trial-level data with per-user ES columns).
"""

import pandas as pd
import pickle


def load_data(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)


def clip_extreme_rt(df, lower=0.3, upper=12.0):
    """Clip response times to plausible range per game type."""
    for _, group in df.groupby('gameName'):
        df.loc[group.index, 'responseTime'] = group['responseTime'].clip(
            lower=lower, upper=upper)
    return df


def compute_es_per_game(df):
    """Compute Efficacy Score (SR / RT) aggregated per gameID."""
    es = df.groupby('gameID').agg(
        mRT=('responseTime', 'mean'),
        mSR=('taskSuccessRate', 'mean')
    ).reset_index()
    es['ES'] = es['mSR'] / es['mRT']
    return df.merge(es[['gameID', 'ES']], on='gameID', how='left')


def add_es_columns(df):
    """Pivot first-session ES per game type into user-level columns."""
    df['startDate'] = pd.to_datetime(df['startDate'])
    df = df.sort_values(by=['userID', 'gameName', 'startDate'])
    df_first = df.drop_duplicates(subset=['userID', 'gameName'], keep='first')

    pivot = df_first.pivot_table(index='userID', columns='gameName',
                                 values='ES', aggfunc='first')
    pivot.columns = [f'ES_{col}' for col in pivot.columns]
    pivot = pivot.reset_index()

    return df.merge(pivot, on='userID', how='left')


def main(input_file='merged_data.pkl', output_dir='.'):
    import os
    os.makedirs(output_dir, exist_ok=True)

    df = load_data(input_file)
    df = clip_extreme_rt(df)
    df = compute_es_per_game(df)
    df = add_es_columns(df)

    pkl_path = os.path.join(output_dir, 'merged_df_ES.pkl')
    csv_path = os.path.join(output_dir, 'merged_df_ES.csv')

    with open(pkl_path, 'wb') as f:
        pickle.dump(df, f)
    df.to_csv(csv_path, index=False)

    print(f"Saved merged_df_ES.pkl and .csv — {df['userID'].nunique()} unique users")
    return df


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Compute Efficacy Scores')
    parser.add_argument('--input', default='merged_data.pkl',
                        help='Path to merged_data.pkl (from preprocessing.py)')
    parser.add_argument('--output-dir', default='.',
                        help='Directory for output files')
    args = parser.parse_args()
    main(args.input, args.output_dir)
