"""
Binary risk-factor definitions based on the 2024 Lancet Commission
on dementia prevention, intervention, and care.

Loads the analysis-ready CSV (merged_df_ES.csv), aggregates to one row
per participant, applies a medical-data availability filter, and creates
13 binary risk-factor columns plus an air-pollution variable.

Risk factors and their definitions
-----------------------------------
education_risk          : Years of education <= 12
hearing_loss_risk       : Self-reported hearing loss without hearing-aid use
high_cholesterol_risk   : Self-reported high cholesterol diagnosis
depression_risk         : Current diagnosis or >= 5 years of depression
diabetes_risk           : Self-reported diabetes diagnosis
smoking_risk            : Active/Heavy Exposure (AHE) —
                          current smoker >= 0.5 packs/day, OR
                          former smoker with >= 10 pack-years
hypertension_risk       : Self-reported hypertension diagnosis
obesity_risk            : BMI >= 30
excessive_alcohol_risk  : >= 14 drinks/week (men) or >= 7 drinks/week (women),
                          following Israeli Ministry of Health guidelines
physical_inactivity_risk: No aerobic exercise OR zero exercise sessions/week
social_isolation_risk   : Lives alone, perceived social decline, or no social meetings
sleep_disorders_risk    : Self-reported insomnia, or sleep < 6 h or > 9 h per night
anxiety_risk            : Current diagnosis or >= 5 years of anxiety
"""

import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')


class RiskFactorDefinitions:
    """Create binary Lancet-Commission risk factors from questionnaire data."""

    def __init__(self):
        self.risk_factor_mapping = {
            'education': {
                'columns': ['demog_edu'],
                'description': 'Less education (5 % risk reduction if improved)',
                'create_binary': self._education,
            },
            'hearing_loss': {
                'columns': ['hear_loss', 'hear_aid'],
                'description': 'Hearing loss (7 % risk reduction if prevented)',
                'create_binary': self._hearing_loss,
            },
            'high_cholesterol': {
                'columns': ['chol_high', 'chol_medFood'],
                'description': 'High LDL cholesterol (7 % risk reduction if managed)',
                'create_binary': self._cholesterol,
            },
            'depression': {
                'columns': ['dep_dia', 'dep_years', 'dep_age'],
                'description': 'Depression (3 % risk reduction if treated)',
                'create_binary': self._depression,
            },
            'diabetes': {
                'columns': ['diab_have', 'diab_dis'],
                'description': 'Diabetes (2 % risk reduction if prevented)',
                'create_binary': self._diabetes,
            },
            'smoking': {
                'columns': ['smoke_now', 'smoke_past', 'smoke_years',
                            'smoke_perDay_f', 'smoke_perDay_m_v2',
                            'smoke_perDay_m_ver2', 'smoke_perDay_past'],
                'description': 'Smoking AHE (2 % risk reduction if eliminated)',
                'create_binary': self._smoking_ahe,
            },
            'hypertension': {
                'columns': ['card_hypertension', 'card_med'],
                'description': 'Hypertension (2 % risk reduction if controlled)',
                'create_binary': self._hypertension,
            },
            'obesity': {
                'columns': ['BMI'],
                'description': 'Obesity (1 % risk reduction if prevented)',
                'create_binary': lambda df: self._simple(df, 'BMI',
                                                         lambda x: (x >= 30).astype(float)),
            },
            'excessive_alcohol': {
                'columns': ['alch_week', 'alch_perWeek', 'gender'],
                'description': 'Excessive alcohol (1 % risk reduction if reduced)',
                'create_binary': self._alcohol,
            },
            'physical_inactivity': {
                'columns': ['sport_aerobic', 'sport_perWeek'],
                'description': 'Physical inactivity (2 % risk reduction if increased)',
                'create_binary': self._physical_inactivity,
            },
            'social_isolation': {
                'columns': ['soc_alone', 'soc_desDecline', 'soc_meetings'],
                'description': 'Social isolation (5 % risk reduction if addressed)',
                'create_binary': self._social_isolation,
            },
            'sleep_disorders': {
                'columns': ['sleep_insom', 'sleep_hours'],
                'description': 'Sleep disorders',
                'create_binary': self._sleep_disorders,
            },
            'anxiety': {
                'columns': ['anx_dia', 'anx_years', 'anx_age'],
                'description': 'Anxiety disorders',
                'create_binary': self._anxiety,
            },
        }
        self.es_columns = ['ES_Places', 'ES_People', 'ES_Events', 'ES_Lexical']

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    def _has_medical(self, df):
        cols = ['alz_family', 'dem_family', 'card_hypertension', 'diab_have',
                'smoke_now', 'dep_dia', 'anx_dia', 'chol_high', 'hear_loss', 'BMI']
        avail = [c for c in cols if c in df.columns]
        return df[avail].notna().any(axis=1)

    def _simple(self, df, column, func):
        med = self._has_medical(df)
        result = pd.Series(np.nan, index=df.index)
        if column in df.columns:
            valid = med & df[column].notna()
            result.loc[valid] = func(df.loc[valid, column])
        return result

    # ------------------------------------------------------------------
    # Individual risk factors
    # ------------------------------------------------------------------
    def _education(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        if 'demog_edu' in df.columns:
            v = med & df['demog_edu'].notna()
            r.loc[v] = (df.loc[v, 'demog_edu'] <= 12).astype(float)
        return r

    def _hearing_loss(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        if 'hear_loss' in df.columns:
            v = med & (df['hear_loss'].notna() | df.get('hear_aid', pd.Series()).notna())
            loss = df.get('hear_loss', 0).fillna(0)
            aid = df.get('hear_aid', 0).fillna(0)
            r.loc[v] = ((loss == 1) & (aid == 0)).astype(float)
        return r

    def _cholesterol(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        if 'chol_high' in df.columns:
            v = med & df['chol_high'].notna()
            r.loc[v] = (df.loc[v, 'chol_high'] == 1).astype(float)
        return r

    def _depression(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        if 'dep_dia' in df.columns:
            v = med & df['dep_dia'].notna()
            flag = (df.get('dep_dia', 0) == 1).fillna(False)
            if 'dep_years' in df.columns:
                flag = flag | (df.get('dep_years', 0) >= 5).fillna(False)
            r.loc[v] = flag.loc[v].astype(float)
        return r

    def _diabetes(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        cols = [c for c in ['diab_have', 'diab_dis'] if c in df.columns]
        if cols:
            v = med & df[cols].notna().any(axis=1)
            r.loc[v] = (df[cols] == 1).any(axis=1).loc[v].astype(float)
        return r

    def _smoking_ahe(self, df):
        """Active / Heavy Exposure: current >= 0.5 packs/day OR former >= 10 pack-years."""
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        cols = [c for c in ['smoke_now', 'smoke_past'] if c in df.columns]
        if cols:
            v = med & df[cols].notna().any(axis=1)
            current = (df.get('smoke_now', 0) == 1).fillna(False)
            past = (df.get('smoke_past', 0) == 1).fillna(False)
            years = df.get('smoke_years', 0).fillna(0)

            cpd = pd.Series(0.0, index=df.index)
            if 'smoke_perDay_f' in df.columns:
                cpd = cpd.combine_first(df['smoke_perDay_f'].fillna(0) / 20)
            if 'smoke_perDay_m_v2' in df.columns:
                cpd = cpd.combine_first(df['smoke_perDay_m_v2'].fillna(0) / 20)

            past_py = pd.Series(0.0, index=df.index)
            if 'smoke_perDay_past' in df.columns:
                past_py = df['smoke_perDay_past'].fillna(0) / 20 * years

            ahe = (current & (cpd >= 0.5)) | (past & (past_py >= 10))
            r.loc[v] = ahe.loc[v].astype(float)
        return r

    def _hypertension(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        if 'card_hypertension' in df.columns:
            v = med & df['card_hypertension'].notna()
            r.loc[v] = (df.loc[v, 'card_hypertension'] == 1).astype(float)
        return r

    def _alcohol(self, df):
        """Israeli MoH thresholds: >= 14 drinks/week (men), >= 7 (women)."""
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        a_cols = [c for c in ['alch_week', 'alch_perWeek'] if c in df.columns]
        if a_cols and 'gender' in df.columns:
            v = med & df[a_cols].notna().any(axis=1) & df['gender'].notna()
            apw = df.get('alch_perWeek', pd.Series(0.0, index=df.index)).fillna(0)
            if 'alch_week' in df.columns:
                apw = apw.combine_first(df['alch_week'].fillna(0))
            male = df['gender'] == 'MALE'
            female = df['gender'] == 'FEMALE'
            exc = (male & (apw >= 14)) | (female & (apw >= 7)) | (~male & ~female & (apw >= 7))
            r.loc[v] = exc.loc[v].astype(float)
        return r

    def _physical_inactivity(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        cols = [c for c in ['sport_aerobic', 'sport_perWeek'] if c in df.columns]
        if cols:
            v = med & df[cols].notna().any(axis=1)
            no_aer = (df.get('sport_aerobic', 1) == 0).fillna(True)
            no_wk = (df.get('sport_perWeek', 1) == 0).fillna(True)
            r.loc[v] = (no_aer | no_wk).loc[v].astype(float)
        return r

    def _social_isolation(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        cols = [c for c in ['soc_alone', 'soc_desDecline', 'soc_meetings'] if c in df.columns]
        if cols:
            v = med & df[cols].notna().any(axis=1)
            alone = (df.get('soc_alone', 0) == 1).fillna(False)
            decline = (df.get('soc_desDecline', 0) == 1).fillna(False)
            no_meet = (df.get('soc_meetings', 1) == 0).fillna(False)
            r.loc[v] = (alone | decline | no_meet).loc[v].astype(float)
        return r

    def _sleep_disorders(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        cols = [c for c in ['sleep_insom', 'sleep_hours'] if c in df.columns]
        if cols:
            v = med & df[cols].notna().any(axis=1)
            insom = (df.get('sleep_insom', 0) == 1).fillna(False)
            extreme = ((df.get('sleep_hours', 7.5) < 6) |
                       (df.get('sleep_hours', 7.5) > 9)).fillna(False)
            r.loc[v] = (insom | extreme).loc[v].astype(float)
        return r

    def _anxiety(self, df):
        med = self._has_medical(df)
        r = pd.Series(np.nan, index=df.index)
        if 'anx_dia' in df.columns:
            v = med & df['anx_dia'].notna()
            flag = (df.get('anx_dia', 0) == 1).fillna(False)
            if 'anx_years' in df.columns:
                flag = flag | (df.get('anx_years', 0) >= 5).fillna(False)
            r.loc[v] = flag.loc[v].astype(float)
        return r

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def load_data(self, csv_path):
        """Load CSV, aggregate to user level, and filter to medical respondents."""
        self.df = pd.read_csv(csv_path, low_memory=False)
        self.user_df = self.df.groupby('userID').agg({
            **{c: 'first' for c in self.df.columns
               if c not in ['userID'] + self.es_columns},
            **{c: 'mean' for c in self.es_columns if c in self.df.columns},
        }).reset_index()
        med = self._has_medical(self.user_df)
        self.user_df = self.user_df[med].reset_index(drop=True)
        return self.user_df

    def create_risk_categories(self):
        """Create all binary risk columns and return a summary dict."""
        self.risk_categories = {}
        for name, cfg in self.risk_factor_mapping.items():
            self.user_df[f'{name}_risk'] = cfg['create_binary'](self.user_df)
            valid = self.user_df[f'{name}_risk'].notna()
            n_valid = valid.sum()
            n_pos = (self.user_df[f'{name}_risk'] == 1).sum()
            if n_valid:
                self.risk_categories[name] = {
                    'column': f'{name}_risk',
                    'description': cfg['description'],
                    'n_users': int(n_valid),
                    'n_positive': int(n_pos),
                    'percentage': float(n_pos / n_valid * 100),
                }
        return self.risk_categories
