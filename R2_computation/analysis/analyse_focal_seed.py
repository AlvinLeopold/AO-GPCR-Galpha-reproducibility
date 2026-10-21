#!/usr/bin/env python3
"""Compare focal receptor–Gα scores under archived seed 0 and rerun seed 101."""

from pathlib import Path
import json

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
original = pd.read_csv(ROOT / 'inputs/af2_all_model_scores.tsv', sep='\t')
rerun = pd.read_csv(ROOT / 'inputs/af2_focal_seed101_scores.tsv', sep='\t')
focal_prefixes = ('AOL_s00007g309_1__', 'AOL_s00176g42_1__')
original = original.loc[original['job'].str.startswith(focal_prefixes)].copy()
assert len(original) == len(rerun) == 18
assert set(original['seed']) == {0} and set(rerun['seed']) == {101}
combined = pd.concat([original, rerun], ignore_index=True)
combined['receptor'] = combined['job'].str.split('__').str[0].str.replace('_1', '.1', regex=False)
combined['galpha'] = combined['job'].str.split('__').str[1]
assert combined.groupby(['receptor', 'galpha', 'seed']).size().eq(3).all()
assert combined.groupby(['receptor', 'galpha', 'seed'])['model'].nunique().eq(3).all()
combined.sort_values(['receptor', 'seed', 'galpha', 'model']).to_csv(ROOT / 'analysis/af2_focal_two_seed_scores.tsv', sep='\t', index=False)

records = []
for (receptor, seed), d in combined.groupby(['receptor', 'seed']):
    max_score = d.groupby('galpha')['iptm'].max().sort_values(ascending=False)
    mean_score = d.groupby('galpha')['iptm'].mean().sort_values(ascending=False)
    model_winners = []
    for _, sub in d.groupby('model'):
        model_winners.append(sub.loc[sub['iptm'].idxmax(), 'galpha'])
    records.append({
        'receptor': receptor,
        'seed': int(seed),
        'winner_by_max': max_score.index[0],
        'winner_by_mean': mean_score.index[0],
        'max_margin': float(max_score.iloc[0] - max_score.iloc[1]),
        'mean_margin': float(mean_score.iloc[0] - mean_score.iloc[1]),
        'model_1_winner': model_winners[0],
        'model_2_winner': model_winners[1],
        'model_3_winner': model_winners[2],
    })
summary = pd.DataFrame(records).sort_values(['receptor', 'seed'])
summary.to_csv(ROOT / 'analysis/af2_focal_seed_sensitivity.tsv', sep='\t', index=False)

pair_stats = combined.groupby(['receptor', 'galpha', 'seed'])['iptm'].agg(['min', 'mean', 'max']).reset_index()
pair_stats.to_csv(ROOT / 'analysis/af2_focal_seed_pair_ranges.tsv', sep='\t', index=False)
table_s6 = pair_stats.merge(summary, on=['receptor', 'seed'], validate='many_to_one')
(ROOT / 'supplementary').mkdir(exist_ok=True)
table_s6.to_csv(ROOT / 'supplementary/Supplementary_Table_S6_AF2_Seed_Sensitivity.tsv', sep='\t', index=False)
report = {
    'scores_total': len(combined),
    'scores_per_seed': {str(int(k)): int(v) for k,v in combined.groupby('seed').size().items()},
    'max_winner_stable_across_seeds': bool(summary.groupby('receptor')['winner_by_max'].nunique().eq(1).all()),
    'mean_winner_stable_across_seeds': bool(summary.groupby('receptor')['winner_by_mean'].nunique().eq(1).all()),
    'summary': records,
}
(ROOT / 'analysis/af2_focal_seed_summary.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
