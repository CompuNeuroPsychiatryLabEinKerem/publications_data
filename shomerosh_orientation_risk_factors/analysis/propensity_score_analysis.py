"""
Overlap-weighted doubly robust propensity-score analysis.

Propensity score  PS = P(Exposure | Age, Age², Sex, Age × Sex)
is estimated via logistic regression.  Two doubly robust approaches
are implemented:

  Approach A — PS as an additional covariate in OLS:
      Score ~ Exposure + PS + Age + Age² + Sex + Age × Sex

  Approach B — Overlap weights in weighted least squares:
      Weight = (1 − PS)  for exposed participants
      Weight = PS         for unexposed participants
      Score ~ Exposure + Age + Age² + Sex + Age × Sex   [WLS]

The primary outcome is the incremental variance explained (ΔR²)
attributable to the risk factor after adjusting for demographics.

Reference:
  Li, Morgan & Zaslavsky (2018).  Balancing covariates via propensity
  score weighting.  JASA 113(521):390-400.
"""

import os
import sys
import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from risk_factor_definitions import RiskFactorDefinitions

ES_MEASURES = [
    ('ES_Places',               'Spatial (Places)'),
    ('ES_Events',               'Temporal (Events)'),
    ('ES_People',               'Social (People)'),
    ('ES_Lexical',              'Lexical (Language)'),
    ('ES_Spatial_Temporal',     'Spatial-Temporal'),
    ('ES_Combined_Orientation', 'Combined Orientation'),
    ('ES_Overall_Cognitive',    'Overall Cognitive'),
]


def load_and_prepare(csv_path, air_pollution_dir=None):
    """Load data, create risk factors, and add derived columns."""
    analyzer = RiskFactorDefinitions()
    analyzer.load_data(csv_path)
    analyzer.create_risk_categories()
    df = analyzer.user_df.copy()

    df['Sex'] = (df['gender'] == 'MALE').astype(float)
    df['Age_squared'] = df['Age'] ** 2
    df['Age_x_Sex'] = df['Age'] * df['Sex']
    df = df.dropna(subset=['Age', 'gender'])

    df['ES_Spatial_Temporal'] = df[['ES_Places', 'ES_Events']].mean(axis=1)
    df['ES_Combined_Orientation'] = df[['ES_Places', 'ES_Events', 'ES_People']].mean(axis=1)
    df['ES_Overall_Cognitive'] = df[['ES_Places', 'ES_Events', 'ES_People', 'ES_Lexical']].mean(axis=1)

    if air_pollution_dir:
        _add_air_pollution(df, csv_path, air_pollution_dir)

    return df


def _add_air_pollution(df, csv_path, air_dir):
    """Map participant city of residence to Comprehensive Pollution Index."""
    try:
        mapping_path = os.path.join(air_dir, 'city_area_mapping.csv')
        if not os.path.exists(mapping_path):
            return
        raw = pd.read_csv(csv_path, low_memory=False)
        if 'cityName' not in raw.columns:
            return
        city_info = raw.groupby('userID')['cityName'].first().reset_index()
        if 'cityName' in df.columns:
            df.drop(columns=['cityName'], inplace=True)
        df_tmp = df[['userID']].merge(city_info, on='userID', how='left')
        cm = pd.read_csv(mapping_path)
        cpi_map = dict(zip(cm['city'], cm['CPI_mean']))
        df['air_pollution_CPI'] = df_tmp['cityName'].map(cpi_map).values
        df['air_pollution_very_high_risk'] = (df['air_pollution_CPI'] > 2.0).astype(float)
        df.loc[df['air_pollution_CPI'].isna(), 'air_pollution_very_high_risk'] = np.nan
    except Exception:
        pass


def run_propensity_analysis(df, score_col, exposure_col):
    """Run both PS approaches for one score × exposure pair."""
    req = [score_col, exposure_col, 'Age', 'Sex', 'Age_squared', 'Age_x_Sex']
    adf = df[req].dropna().copy()
    if len(adf) < 50 or adf[exposure_col].nunique() < 2:
        return None
    n_exp = (adf[exposure_col] == 1).sum()
    n_unexp = (adf[exposure_col] == 0).sum()
    if n_exp < 10 or n_unexp < 10:
        return None

    adf = adf.rename(columns={score_col: 'Score', exposure_col: 'Exposure'})

    try:
        # Step 1 — Propensity score
        ps_model = smf.logit(
            'Exposure ~ Age + Age_squared + Sex + Age_x_Sex', data=adf
        ).fit(disp=0, method='bfgs', maxiter=1000)
        adf['PS'] = ps_model.predict()

        # Null (demographics only)
        null_formula = 'Score ~ Age + Age_squared + Sex + Age_x_Sex'
        m_null = smf.ols(null_formula, data=adf).fit()

        # Approach A — PS as covariate
        m_a = smf.ols(
            'Score ~ Exposure + PS + Age + Age_squared + Sex + Age_x_Sex', data=adf
        ).fit()

        # Approach B — Overlap weights
        adf['ow'] = np.where(adf['Exposure'] == 1, 1 - adf['PS'], adf['PS'])
        m_null_w = smf.wls(null_formula, data=adf, weights=adf['ow']).fit()
        m_b = smf.wls(
            'Score ~ Exposure + Age + Age_squared + Sex + Age_x_Sex',
            data=adf, weights=adf['ow']
        ).fit()

        return {
            'score': score_col,
            'exposure': exposure_col,
            'n': len(adf), 'n_exposed': int(n_exp), 'n_unexposed': int(n_unexp),
            'r2_null': m_null.rsquared,
            # Approach A
            'dr2_A': (m_a.rsquared - m_null.rsquared) * 100,
            'beta_A': m_a.params.get('Exposure', np.nan),
            'p_A': m_a.pvalues.get('Exposure', np.nan),
            'sig_A': m_a.pvalues.get('Exposure', 1) < 0.05,
            # Approach B
            'dr2_B': (m_b.rsquared - m_null_w.rsquared) * 100,
            'beta_B': m_b.params.get('Exposure', np.nan),
            'p_B': m_b.pvalues.get('Exposure', np.nan),
            'sig_B': m_b.pvalues.get('Exposure', 1) < 0.05,
        }
    except Exception:
        return None


def analyse_all(df):
    """Run both approaches for every risk factor × cognitive measure."""
    rfs = sorted(c for c in df.columns if c.endswith('_risk'))
    rows = []
    for es_col, es_name in ES_MEASURES:
        if es_col not in df.columns:
            continue
        for rf in rfs:
            r = run_propensity_analysis(df, es_col, rf)
            if r:
                r['display_name'] = es_name
                rows.append(r)
    return pd.DataFrame(rows)


def main(csv_path, air_pollution_dir=None, output_dir='.'):
    os.makedirs(output_dir, exist_ok=True)
    df = load_and_prepare(csv_path, air_pollution_dir)
    results = analyse_all(df)
    results.to_csv(os.path.join(output_dir, 'propensity_results.csv'), index=False)
    print(f"Saved propensity_results.csv — {len(results)} tests")
    return results, df


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Overlap-weighted doubly robust propensity-score analysis')
    parser.add_argument('--data', default='../preprocessing/merged_df_ES.csv',
                        help='Path to merged_df_ES.csv (from efficacy_scoring.py)')
    parser.add_argument('--air-pollution-dir', default=None,
                        help='Directory with city_area_mapping.csv (optional)')
    parser.add_argument('--output-dir', default='results',
                        help='Directory for output files')
    args = parser.parse_args()
    main(args.data, args.air_pollution_dir, args.output_dir)
