"""
Analysis pipeline for "Mid-Infrared Colours Trace Formation Channel in the
COUGS-DESI Catalogue of Polar-Structure Galaxies" (Dufera 2026).

This script reproduces the full analysis: cross-matching the COUGS-DESI
catalogue against AllWISE, computing WISE colours, and performing every
statistical test reported in the Letter.

Data sources (both public, no authentication required):
    COUGS-DESI: VizieR J/A+A/710/A145 (Bahr-Kalus & Mosenkov 2026)
    AllWISE:    VizieR II/328/allwise (Cutri et al. 2013)
"""
import numpy as np
import pandas as pd
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.stats import kruskal, mannwhitneyu, spearmanr, fisher_exact, chi2_contingency

STERN_W1W2_THRESHOLD = 0.8
STERN_W2_MAG_LIMIT = 15.05
QUIESCENT_W2W3_THRESHOLD = 2.0
STARFORMING_W2W3_THRESHOLD = 3.5
WISE_MATCH_RADIUS_ARCSEC = 6


def load_cougs_desi():
    """Pull the real COUGS-DESI catalogue and remove placeholder redshifts."""
    Vizier.ROW_LIMIT = -1
    result = Vizier.get_catalogs("J/A+A/710/A145/cougs")
    df = result[0].to_pandas()
    clean = df[df['z'] != -1.0].dropna(subset=['z']).copy()
    return clean


def crossmatch_wise(cougs_df):
    """
    Cross-match every COUGS-DESI galaxy against AllWISE, retaining the
    nearest match within WISE_MATCH_RADIUS_ARCSEC.
    """
    Vizier.ROW_LIMIT = 5
    Vizier.columns = ['RAJ2000', 'DEJ2000', 'W1mag', 'W2mag', 'W3mag', 'W4mag']

    results = []
    for _, target in cougs_df.iterrows():
        target_coord = SkyCoord(ra=target['RAJ2000'] * u.deg, dec=target['DEJ2000'] * u.deg)
        try:
            result = Vizier.query_region(
                target_coord, radius=WISE_MATCH_RADIUS_ARCSEC * u.arcsec, catalog='II/328/allwise')
            if len(result) > 0 and len(result[0]) > 0:
                row = result[0][0]
                results.append({
                    'Name': target['Name'], 'TYPE1': target['TYPE1'], 'z': target['z'],
                    'W1mag': row['W1mag'], 'W2mag': row['W2mag'],
                    'W3mag': row['W3mag'], 'W4mag': row['W4mag'], 'matched': True,
                })
            else:
                results.append({'Name': target['Name'], 'TYPE1': target['TYPE1'], 'z': target['z'],
                                 'W1mag': None, 'W2mag': None, 'W3mag': None, 'W4mag': None, 'matched': False})
        except Exception:
            results.append({'Name': target['Name'], 'TYPE1': target['TYPE1'], 'z': target['z'],
                             'W1mag': None, 'W2mag': None, 'W3mag': None, 'W4mag': None, 'matched': False})

    wise_df = pd.DataFrame(results)
    for col in ['W1mag', 'W2mag', 'W3mag', 'W4mag']:
        wise_df[col] = pd.to_numeric(wise_df[col], errors='coerce').astype('float64')

    wise_df['W1_W2'] = wise_df['W1mag'] - wise_df['W2mag']
    wise_df['W2_W3'] = wise_df['W2mag'] - wise_df['W3mag']
    wise_df['W1_W4'] = wise_df['W1mag'] - wise_df['W4mag']
    return wise_df


def classify_wise(row):
    """
    Real, literature-correct three-zone WISE classification (Jarrett et al.
    2011, 2017), with the Stern et al. (2012) AGN criterion including its
    magnitude qualifier.
    """
    if pd.isna(row['W2_W3']):
        return None
    if row['W1_W2'] > STERN_W1W2_THRESHOLD and row['W2mag'] < STERN_W2_MAG_LIMIT:
        return 'AGN'
    elif row['W2_W3'] < QUIESCENT_W2W3_THRESHOLD:
        return 'Quiescent'
    elif row['W2_W3'] <= STARFORMING_W2W3_THRESHOLD:
        return 'Intermediate'
    else:
        return 'Star-forming'


def bootstrap_median_diff(group_a, group_b, n_boot=5000, seed=42):
    """Bootstrap confidence interval on the difference of two group medians."""
    rng = np.random.default_rng(seed)
    diffs = np.array([
        np.median(rng.choice(group_b, size=len(group_b), replace=True))
        - np.median(rng.choice(group_a, size=len(group_a), replace=True))
        for _ in range(n_boot)
    ])
    return np.median(group_b) - np.median(group_a), np.percentile(diffs, 2.5), np.percentile(diffs, 97.5)


def run_full_analysis():
    print("Loading COUGS-DESI...")
    cougs = load_cougs_desi()
    print(f"  {len(cougs)} usable galaxies")

    print("Cross-matching against AllWISE (this takes a while for the full sample)...")
    wise_df = crossmatch_wise(cougs)
    matched = wise_df[wise_df['matched']].dropna(subset=['W1_W2']).copy()
    print(f"  Matched: {matched['matched'].sum()} / {len(wise_df)}")

    matched_w23 = matched.dropna(subset=['W2_W3']).copy()
    print(f"  Usable for W2-W3 (W3-detected): {len(matched_w23)}")

    # --- Section 4.1: overall subtype differences ---
    print("\n--- WISE colours by subtype ---")
    print(matched_w23.groupby('TYPE1')[['W1_W2', 'W2_W3']].median())

    groups_w1w2 = [g['W1_W2'].values for _, g in matched.groupby('TYPE1')]
    groups_w2w3 = [g['W2_W3'].values for _, g in matched_w23.groupby('TYPE1')]
    h1, p1 = kruskal(*groups_w1w2)
    h2, p2 = kruskal(*groups_w2w3)
    print(f"Kruskal-Wallis W1-W2: H={h1:.2f}, p={p1:.4f}")
    print(f"Kruskal-Wallis W2-W3: H={h2:.2f}, p={p2:.4f}")

    rho1, prho1 = spearmanr(matched['z'], matched['W1_W2'])
    rho2, prho2 = spearmanr(matched_w23['z'], matched_w23['W2_W3'])
    print(f"Redshift correlation, W1-W2: rho={rho1:.3f}, p={prho1:.4f}")
    print(f"Redshift correlation, W2-W3: rho={rho2:.3f}, p={prho2:.4f}")

    # --- Section 4.2: PR vs PTS ---
    print("\n--- PR vs PTS (core prediction) ---")
    pr = matched_w23.loc[matched_w23['TYPE1'] == 'PR', 'W2_W3'].values
    pts = matched_w23.loc[matched_w23['TYPE1'] == 'PTS', 'W2_W3'].values
    u_stat, p_prpts = mannwhitneyu(pr, pts, alternative='two-sided')
    diff, lo, hi = bootstrap_median_diff(pr, pts)
    print(f"PR median: {np.median(pr):.3f}, PTS median: {np.median(pts):.3f}, p={p_prpts:.6f}")
    print(f"Bootstrap difference: {diff:.3f}, 95% CI: [{lo:.3f}, {hi:.3f}]")

    # --- Section 4.3: PB as outlier, correct 3-zone classification ---
    print("\n--- PB subtype analysis ---")
    matched_w23['wise_class'] = matched_w23.apply(classify_wise, axis=1)

    pb_mask = matched_w23['TYPE1'] == 'PB'
    pb_q = (matched_w23.loc[pb_mask, 'wise_class'] == 'Quiescent').sum()
    other_q = (matched_w23.loc[~pb_mask, 'wise_class'] == 'Quiescent').sum()
    n_pb, n_other = pb_mask.sum(), (~pb_mask).sum()

    print(f"PB quiescent fraction: {pb_q}/{n_pb} = {pb_q/n_pb*100:.1f}%")
    print(f"Other subtypes quiescent fraction: {other_q}/{n_other} = {other_q/n_other*100:.1f}%")
    _, p_q = fisher_exact([[pb_q, n_pb - pb_q], [other_q, n_other - other_q]])
    print(f"Fisher's exact test (quiescent): p={p_q:.6f}")

    pb_sf = (matched_w23.loc[pb_mask, 'wise_class'] == 'Star-forming').sum()
    other_sf = (matched_w23.loc[~pb_mask, 'wise_class'] == 'Star-forming').sum()
    _, p_sf = fisher_exact([[pb_sf, n_pb - pb_sf], [other_sf, n_other - other_sf]])
    print(f"PB star-forming fraction: {pb_sf}/{n_pb} = {pb_sf/n_pb*100:.1f}%")
    print(f"Other subtypes star-forming fraction: {other_sf}/{n_other} = {other_sf/n_other*100:.1f}%")
    print(f"Fisher's exact test (star-forming): p={p_sf:.6f}")

    # --- Section 4.4: AGN fraction ---
    print("\n--- AGN fraction ---")
    n_agn = (matched_w23['wise_class'] == 'AGN').sum()
    print(f"AGN-classified: {n_agn} / {len(matched_w23)}")
    contingency = pd.crosstab(matched_w23['TYPE1'], matched_w23['wise_class'] == 'AGN')
    chi2, p_agn, _, _ = chi2_contingency(contingency)
    print(f"Chi-squared test, AGN fraction across subtypes: chi2={chi2:.2f}, p={p_agn:.4f}")

    return matched_w23


if __name__ == '__main__':
    run_full_analysis()
