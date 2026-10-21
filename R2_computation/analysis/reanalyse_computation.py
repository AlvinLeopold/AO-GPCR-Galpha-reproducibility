#!/usr/bin/env python3
"""Auditable sensitivity analyses for IJMS-4611842 computational revision.

Uses only archived source tables and the per-model score JSON extraction.
No prediction is rerun and no new biological replicate is inferred.
"""

from pathlib import Path
import csv
import json
import math
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "inputs"
OUT = ROOT / "analysis"
GEO = IN / "atcc_gpcr_analysis" / "results"
OUT.mkdir(exist_ok=True)


def save(df, name):
    df.to_csv(OUT / name, sep="\t", index=False, na_rep="NA")


def main():
    selection = pd.read_csv(IN / "selection_criterion.tsv", sep="\t")
    panel = pd.read_csv(IN / "Figure2_SourceData_PanelC.tsv", sep="\t")
    effects = pd.read_csv(GEO / "effects" / "candidate_effects.tsv", sep="\t")
    evidence = pd.read_csv(GEO / "integration" / "candidate_evidence_by_gse.tsv", sep="\t")
    inventory = pd.read_csv(IN / "tree_annotation_metadata.tsv", sep="\t")
    af2 = pd.read_csv(IN / "af2_all_model_scores.tsv", sep="\t")

    assert len(selection) == 16 and len(panel) == 13
    assert selection.candidate_id.nunique() == 16
    assert af2.job.nunique() == 267 and len(af2) == 801
    assert set(af2.groupby("job").size()) == {3}
    assert set(af2.seed) == {0}

    # Reconstruct the submitted 13-candidate heuristic directly from the
    # contrast-level processed expression table: median across time points for
    # each regulator, then maximum absolute regulator effect in each module.
    matrix_to_regulator = {
        "GSE210886_Acy": "AoAcy",
        "GSE210886_PKA": "AoPkaC1",
        "GSE306156_TPM": "AoOpy2",
        "GSE322548_TPM": "AoHog1",
    }
    x = effects[
        effects.candidate_id.isin(selection.candidate_id)
        & effects.contrast_type.eq("genotype_effect")
        & effects.matrix_id.isin(matrix_to_regulator)
    ].copy()
    x["regulator"] = x.matrix_id.map(matrix_to_regulator)
    x["module"] = np.where(x.regulator.isin(["AoAcy", "AoPkaC1"]), "cAMP_PKA", "MAPK")
    regulators = (
        x.groupby(["candidate_id", "module", "regulator"], as_index=False)
        .agg(median_log2fc=("log2fc", "median"), n_contrasts=("log2fc", "size"),
             supported_contrasts=("supported", "sum"))
    )
    regulators["abs_median"] = regulators.median_log2fc.abs()
    save(regulators, "transcriptome_regulator_effects.tsv")

    pivot = regulators.pivot(index="candidate_id", columns="regulator", values="median_log2fc")
    assert not pivot.loc[panel.candidate_id].isna().any().any()
    base = pd.DataFrame(index=panel.candidate_id)
    base["MAPK_max_abs"] = pivot.loc[base.index, ["AoHog1", "AoOpy2"]].abs().max(axis=1)
    base["cAMP_max_abs"] = pivot.loc[base.index, ["AoAcy", "AoPkaC1"]].abs().max(axis=1)
    base["score_reconstructed"] = base.MAPK_max_abs + base.cAMP_max_abs
    base["score_submitted"] = panel.set_index("candidate_id").loc[base.index, "sum"]
    base["difference"] = base.score_reconstructed - base.score_submitted
    # Submitted table rounded each component to three decimal places.
    assert base.difference.abs().max() < 0.0011, base.sort_values("difference").tail()

    base["rank_submitted"] = base.score_submitted.rank(ascending=False, method="min").astype(int)
    base["score_module_median_abs"] = (
        pivot.loc[base.index, ["AoHog1", "AoOpy2"]].abs().median(axis=1)
        + pivot.loc[base.index, ["AoAcy", "AoPkaC1"]].abs().median(axis=1)
    )
    base["rank_module_median_abs"] = base.score_module_median_abs.rank(ascending=False, method="min").astype(int)
    for mapk in ("AoHog1", "AoOpy2"):
        for camp in ("AoAcy", "AoPkaC1"):
            col = f"rank_{mapk}_{camp}"
            score = pivot.loc[base.index, mapk].abs() + pivot.loc[base.index, camp].abs()
            base[col] = score.rank(ascending=False, method="min").astype(int)
    base = base.reset_index().rename(columns={"index": "candidate_id"})
    base = base.sort_values("rank_submitted")
    save(base, "transcriptome_rank_sensitivity.tsv")

    # Sensitivity of the 'WT responsiveness' eligibility component to the two
    # flagged series. The MAPK and cAMP contrasts are from other series.
    e = evidence[evidence.candidate_id.isin(inventory.gpcr)].copy()
    e["support"] = pd.to_numeric(e.wt_induction_supported, errors="coerce").fillna(0) > 0
    rows = []
    for label, removed in [
        ("all_11_series", set()),
        ("without_GSE303416", {"GSE303416"}),
        ("without_GSE253759", {"GSE253759"}),
        ("without_both", {"GSE303416", "GSE253759"}),
    ]:
        sub = e[~e.gse.isin(removed)]
        counts = sub.groupby("candidate_id").support.sum()
        for cid in selection.candidate_id:
            rows.append({"scenario": label, "candidate_id": cid,
                         "wt_supported_series": int(counts.get(cid, 0)),
                         "retains_wt_support": bool(counts.get(cid, 0) > 0),
                         "uncharacterized": cid in set(panel.candidate_id)})
    drop = pd.DataFrame(rows)
    save(drop, "transcriptome_exclusion_sensitivity.tsv")
    leave_one_rows = []
    for excluded in sorted(e.gse.unique()):
        counts = e.loc[~e.gse.eq(excluded)].groupby('candidate_id').support.sum()
        for cid in selection.candidate_id:
            leave_one_rows.append({
                'excluded_series': excluded,
                'candidate_id': cid,
                'wt_supported_series': int(counts.get(cid, 0)),
                'retains_wt_support': bool(counts.get(cid, 0) > 0),
            })
    leave_one = pd.DataFrame(leave_one_rows)
    save(leave_one, 'transcriptome_leave_one_series_out.tsv')
    global_support = {}
    for label, removed in [
        ("all_11_series", set()),
        ("without_GSE303416", {"GSE303416"}),
        ("without_GSE253759", {"GSE253759"}),
        ("without_both", {"GSE303416", "GSE253759"}),
    ]:
        per_receptor = (
            e[~e.gse.isin(removed)]
            .groupby("candidate_id").support.sum()
        )
        contrast_counts = (
            evidence[evidence.candidate_id.isin(inventory.gpcr) & ~evidence.gse.isin(removed)]
            .groupby("candidate_id").wt_induction_supported.sum()
        )
        global_support[label] = {
            "mapped_receptors": int(len(per_receptor)),
            "at_least_one_supported_series": int((per_receptor > 0).sum()),
            "at_least_three_supported_series": int((per_receptor >= 3).sum()),
            "at_least_three_supported_contrasts": int((contrast_counts >= 3).sum()),
        }

    # Re-evaluate the near-zero Kuo expression example using the rounded
    # values reported in the submitted manuscript, and label it approximate.
    pseudo = pd.DataFrame({
        "pseudocount_tpm": [0.01, 0.1, 0.5, 1.0],
        "baseline_tpm_rounded": 0.02,
        "hpe10_tpm_rounded": 1.55,
    })
    pseudo["approx_log2fc"] = np.log2((1.55 + pseudo.pseudocount_tpm) / (0.02 + pseudo.pseudocount_tpm))
    save(pseudo, "low_expression_pseudocount_sensitivity.tsv")

    # AF2: compare partner calls based on maximum ipTM, mean ipTM, and the
    # three distinct model weights. These are NOT three independent seeds.
    parsed = af2.job.str.extract(r"^(AOL_.+)_([0-9]+)__AOGpa([123])$")
    assert parsed.notna().all().all()
    af2["gpcr"] = parsed[0] + "." + parsed[1]
    af2["galpha"] = "AOGpa" + parsed[2]
    assert af2.gpcr.nunique() == 89
    pair = (
        af2.groupby(["gpcr", "galpha"], as_index=False)
        .agg(iptm_min=("iptm", "min"), iptm_median=("iptm", "median"),
             iptm_mean=("iptm", "mean"), iptm_max=("iptm", "max"),
             iptm_sd=("iptm", "std"), n_models=("iptm", "size"))
    )
    save(pair, "af2_pair_all_model_summary.tsv")
    save(af2[["gpcr", "galpha", "model", "seed", "rank", "iptm", "ptm", "max_pae", "file"]],
         "af2_all_model_scores_clean.tsv")

    results = []
    for gpcr, group in pair.groupby("gpcr"):
        by_model = af2[af2.gpcr.eq(gpcr)].pivot(index="model", columns="galpha", values="iptm")
        assert by_model.shape == (3, 3)
        max_rank = group.sort_values(["iptm_max", "galpha"], ascending=[False, True])
        mean_rank = group.sort_values(["iptm_mean", "galpha"], ascending=[False, True])
        winner_max = max_rank.iloc[0].galpha
        winner_mean = mean_rank.iloc[0].galpha
        model_winners = by_model.idxmax(axis=1).tolist()
        out = {"gpcr": gpcr, "winner_max": winner_max, "winner_mean": winner_mean,
               "max_margin": float(max_rank.iloc[0].iptm_max - max_rank.iloc[1].iptm_max),
               "mean_margin": float(mean_rank.iloc[0].iptm_mean - mean_rank.iloc[1].iptm_mean),
               "max_mean_agree": winner_max == winner_mean,
               "same_winner_all_3_model_weights": len(set(model_winners)) == 1,
               "max_winner_supported_by_model_weights": model_winners.count(winner_max),
               "model_1_winner": model_winners[0], "model_2_winner": model_winners[1],
               "model_3_winner": model_winners[2]}
        results.append(out)
    receptors = pd.DataFrame(results)
    save(receptors, "af2_receptor_rank_robustness.tsv")

    submitted = pd.read_csv(IN / "submitted_coupling_selectivity.tsv", sep="\t")
    check = receptors.merge(submitted[["gpcr", "best_galpha", "call"]], on="gpcr", validate="one_to_one")
    assert (check.winner_max == check.best_galpha).all()
    focal = ["AOL_s00007g309.1", "AOL_s00176g42.1", "AOL_s00054g150.1", "AOL_s00076g450.1"]
    focal_pairs = pair[pair.gpcr.isin(focal)].merge(receptors, on="gpcr")
    save(focal_pairs, "af2_focal_pair_sensitivity.tsv")

    summary = {
        "transcriptome": {
            "submitted_selection_n": int(len(selection)),
            "uncharacterized_ranked_n": int(len(panel)),
            "max_reconstruction_error_from_rounded_submitted_score": float(base.difference.abs().max()),
            "drop_scenarios": drop.groupby("scenario").retains_wt_support.sum().astype(int).to_dict(),
            "leave_one_series_out_min_filter_positive_retaining_wt_support": int(leave_one.groupby('excluded_series').retains_wt_support.sum().min()),
            "global_support": global_support,
            "original_top4": base.sort_values("rank_submitted").candidate_id.head(4).tolist(),
            "median_abs_top4": base.sort_values("rank_module_median_abs").candidate_id.head(4).tolist(),
            "focal_alternative_ranks": base[base.candidate_id.isin(focal)][["candidate_id", "rank_submitted", "rank_module_median_abs", "rank_AoHog1_AoAcy", "rank_AoHog1_AoPkaC1", "rank_AoOpy2_AoAcy", "rank_AoOpy2_AoPkaC1"]].to_dict("records"),
        },
        "af2": {
            "pairs": int(len(pair)), "models": int(len(af2)), "seeds": sorted(af2.seed.unique().tolist()),
            "winner_max_counts": receptors.winner_max.value_counts().to_dict(),
            "winner_mean_counts": receptors.winner_mean.value_counts().to_dict(),
            "max_mean_agreement_n": int(receptors.max_mean_agree.sum()),
            "same_winner_all_3_model_weights_n": int(receptors.same_winner_all_3_model_weights.sum()),
            "max_margin_under_0_02_n": int((receptors.max_margin < 0.02).sum()),
            "max_margin_under_0_05_n": int((receptors.max_margin < 0.05).sum()),
            "focal": receptors[receptors.gpcr.isin(focal)].to_dict("records"),
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
