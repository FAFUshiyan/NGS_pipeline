#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ---------- Variant effect sets (Sequence Ontology names from SnpEff) ----------
LOF_EFFECTS = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "exon_loss_variant",
    "transcript_ablation",
}
CANONICAL_SPLICE = {"splice_acceptor_variant", "splice_donor_variant"}
SPLICE_RELATED = {"splice_acceptor_variant", "splice_donor_variant", "splice_region_variant"}
MISSENSE_RELATED = {"missense_variant", "protein_altering_variant", "inframe_insertion", "inframe_deletion"}
BENIGN_LIKE = {
    "synonymous_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intergenic_region",
    "intron_variant",
    "UTR_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_variant",
}


# ---------- Helpers ----------
def sniff_sep(path: str) -> str:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        sample = f.read(8192)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=["\t", ",", ";"])
        return dialect.delimiter
    except Exception:
        return "\t"


def as_float(x):
    if x is None:
        return np.nan
    x = str(x).strip()
    if x in {"", ".", "NA", "nan", "NaN"}:
        return np.nan
    try:
        return float(x)
    except Exception:
        return np.nan


def as_int(x):
    if x is None:
        return np.nan
    x = str(x).strip()
    if x in {"", ".", "NA", "nan", "NaN"}:
        return np.nan
    try:
        return int(float(x))
    except Exception:
        return np.nan


def split_effects(effect_str: str) -> set:
    """
    SnpEff sometimes joins multiple effects like:
      splice_region_variant&intron_variant
    We split on '&' and return a set.
    """
    if effect_str in {None, ".", ""}:
        return set()
    return set(str(effect_str).split("&"))


def parse_gt(gt: str) -> str:
    """
    Return genotype class: HOM_REF / HET / HOM_ALT / OTHER / MISSING
    Accept phased '0|1' too.
    """
    if gt in {None, ".", "", "./.", ".|."}:
        return "MISSING"
    gt = str(gt)
    sep = "/" if "/" in gt else ("|" if "|" in gt else None)
    if not sep:
        return "OTHER"
    a, b = gt.split(sep)
    if a == "." or b == ".":
        return "MISSING"
    if a == "0" and b == "0":
        return "HOM_REF"
    if a == b and a != "0":
        return "HOM_ALT"
    if (a == "0" and b != "0") or (a != "0" and b == "0"):
        return "HET"
    # e.g., 1/2 (multi-allelic) -> treat as HET
    return "HET"


def parse_ad(ad: str):
    """
    Parse AD like 'ref,alt' (or 'ref,alt1,alt2').
    Return (ref_depth, alt_depth_max, vaf_max).
    """
    if ad in {None, ".", ""}:
        return (np.nan, np.nan, np.nan)
    parts = str(ad).split(",")
    vals = [as_float(p) for p in parts]
    if len(vals) < 2 or np.isnan(vals[0]):
        return (np.nan, np.nan, np.nan)
    ref = vals[0]
    alts = [v for v in vals[1:] if not np.isnan(v)]
    if not alts:
        return (ref, np.nan, np.nan)
    alt_max = max(alts)
    denom = ref + alt_max
    vaf = alt_max / denom if denom and not np.isnan(denom) else np.nan
    return (ref, alt_max, vaf)


def compute_af_max(df: pd.DataFrame) -> pd.Series:
    """
    AF_max across all frequency columns available:
    - gnomAD_AF_*
    - WBBC_AF*
    Extend here if you add other population resources.
    """
    af_cols = [c for c in df.columns if c.startswith("gnomAD_AF_") or c.startswith("WBBC_AF")]
    if not af_cols:
        return pd.Series([np.nan] * len(df), index=df.index)

    af = df[af_cols].apply(lambda col: col.map(as_float))
    return af.max(axis=1, skipna=True)


def tier_and_reason(row, args):
    """
    Decide tier (A/B/C/FILTERED) and reasons.
    """
    reasons = []

    gene = row.get(args.col_gene, ".")
    effect_raw = row.get(args.col_effect, ".")
    impact = str(row.get(args.col_impact, "."))
    biotype = str(row.get(args.col_biotype, "."))

    effects = split_effects(effect_raw)
    is_lof = (effects & LOF_EFFECTS) or (impact == "HIGH")
    is_canonical_splice = bool(effects & CANONICAL_SPLICE)
    is_splice_related = bool(effects & SPLICE_RELATED)
    is_missense = bool(effects & MISSENSE_RELATED) or (impact == "MODERATE")

    # scores
    spliceai = as_float(row.get(args.col_spliceai, "."))
    revel = as_float(row.get(args.col_revel, "."))
    cadd = as_float(row.get(args.col_cadd, "."))
    am_score = as_float(row.get(args.col_am_score, "."))
    am_class = str(row.get(args.col_am_class, ".")).lower()

    splice_strong = (not np.isnan(spliceai)) and spliceai >= args.spliceai_strong
    splice_support = (not np.isnan(spliceai)) and spliceai >= args.spliceai_support

    revel_strong = (not np.isnan(revel)) and revel >= args.revel_strong
    revel_support = (not np.isnan(revel)) and revel >= args.revel_support

    cadd_strong = (not np.isnan(cadd)) and cadd >= args.cadd_strong
    cadd_support = (not np.isnan(cadd)) and cadd >= args.cadd_support

    am_lp = (am_class in {"likely_pathogenic", "pathogenic"})
    am_strong = am_lp or ((not np.isnan(am_score)) and am_score >= args.alphamissense_lp)

    # biotype gate (optional)
    if args.protein_coding_only:
        # allow canonical splice / LoF even if biotype missing,
        # but for missense we usually require protein_coding transcript
        if is_missense and biotype not in {"protein_coding"}:
            return ("FILTERED", ["missense_non_protein_coding"])
        # for splice/lof, do not hard-fail on biotype

    # Functional tiering
    if is_lof:
        reasons.append("LoF_or_HIGH")
        return ("A", reasons)

    if is_canonical_splice:
        reasons.append("canonical_splice")
        return ("A", reasons)

    if is_splice_related and splice_strong:
        reasons.append(f"SpliceAI>={args.spliceai_strong}")
        return ("A", reasons)

    if is_missense and (am_strong or revel_strong or cadd_strong):
        reasons.append("missense_strong_predictors")
        return ("A", reasons)

    if is_splice_related and splice_support:
        reasons.append(f"SpliceAI>={args.spliceai_support}")
        return ("B", reasons)

    if is_missense and (revel_support or cadd_support or am_strong):
        reasons.append("missense_support_predictors")
        return ("B", reasons)

    # drop typical benign-like effects unless user wants broad review
    if not args.keep_modifier:
        if effects and effects.issubset(BENIGN_LIKE) and not splice_support:
            return ("FILTERED", ["modifier_or_benignlike_no_splice_support"])

    # Otherwise keep low tier for manual review
    return ("C", ["kept_for_manual_review"])


def main():
    ap = argparse.ArgumentParser(
        description="Clinical-grade rare variant filtering for SnpSift extractFields TSV (top ANN + genotype)."
    )

    ap.add_argument("--input", required=True, help="Input TSV/CSV from SnpSift extractFields.")
    ap.add_argument("--out-prefix", required=True, help="Output prefix for results.")
    ap.add_argument("--sep", default="auto", help="Delimiter: auto/tab/comma/semicolon.")

    # Which sample genotype to use
    ap.add_argument("--sample", default="0", help="Sample selector for GEN[...] fields. Default 0 => GEN[0].*")

    # Inheritance / frequency defaults
    ap.add_argument("--mode", choices=["AD", "AR", "ANY"], default="AD",
                    help="Inheritance model: affects default MAF & AR compound-het logic.")
    ap.add_argument("--maf", type=float, default=None,
                    help="Max AF_max cutoff. Default: AD=0.001, AR=0.01, ANY=0.001")

    # Genotype QC
    ap.add_argument("--min-dp", type=int, default=10)
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--min-alt-depth", type=int, default=3)
    ap.add_argument("--het-vaf-min", type=float, default=0.2)
    ap.add_argument("--het-vaf-max", type=float, default=0.8)
    ap.add_argument("--hom-vaf-min", type=float, default=0.9)

    # Functional thresholds
    ap.add_argument("--spliceai-support", type=float, default=0.2)
    ap.add_argument("--spliceai-strong", type=float, default=0.5)
    ap.add_argument("--revel-support", type=float, default=0.5)
    ap.add_argument("--revel-strong", type=float, default=0.75)
    ap.add_argument("--cadd-support", type=float, default=20.0)
    ap.add_argument("--cadd-strong", type=float, default=25.0)
    ap.add_argument("--alphamissense-lp", type=float, default=0.564)

    # Behavior switches
    ap.add_argument("--protein-coding-only", action="store_true",
                    help="Focus on protein_coding for missense; still keep LoF/splice when biotype missing.")
    ap.add_argument("--keep-modifier", action="store_true",
                    help="If set, keep MODIFIER effects (e.g. upstream/downstream) as tier C instead of filtering.")
    ap.add_argument("--allow-missing-genotype", action="store_true",
                    help="If set, allow rows without genotype fields; otherwise they are filtered out.")
    ap.add_argument("--emit-comphet", action="store_true",
                    help="In AR mode, output compound-het candidate genes summary.")

    # Column names matching your table
    ap.add_argument("--col-gene", default="ANN[0].GENE")
    ap.add_argument("--col-effect", default="ANN[0].EFFECT")
    ap.add_argument("--col-impact", default="ANN[0].IMPACT")
    ap.add_argument("--col-biotype", default="ANN[0].BIOTYPE")
    ap.add_argument("--col-hgvsc", default="ANN[0].HGVS_C")
    ap.add_argument("--col-hgvsp", default="ANN[0].HGVS_P")
    ap.add_argument("--col-cadd", default="gnomAD_cadd_phred")
    ap.add_argument("--col-revel", default="REVEL_max")
    ap.add_argument("--col-spliceai", default="SpliceAI_DS_max")
    ap.add_argument("--col-am-score", default="am_pathogenicity")
    ap.add_argument("--col-am-class", default="am_class")

    args = ap.parse_args()

    # delimiter
    if args.sep == "auto":
        sep = sniff_sep(args.input)
    else:
        sep = {"tab": "\t", "comma": ",", "semicolon": ";"}.get(args.sep, args.sep)

    # default maf by mode
    if args.maf is None:
        args.maf = 0.01 if args.mode == "AR" else 0.001

    df = pd.read_csv(args.input, sep=sep, dtype=str).fillna(".")

    # Build genotype column names
    sid = str(args.sample)
    gt_col = f"GEN[{sid}].GT"
    dp_col = f"GEN[{sid}].DP"
    ad_col = f"GEN[{sid}].AD"
    gq_col = f"GEN[{sid}].GQ"

    has_gt = gt_col in df.columns

    if not has_gt and not args.allow_missing_genotype:
        raise SystemExit(
            f"Genotype column {gt_col} not found. "
            f"Either re-run extractFields with {gt_col} or use --allow-missing-genotype."
        )

    # Compute AF_max
    df["AF_max"] = compute_af_max(df)

    # Variant key
    df["VARIANT_KEY"] = df["CHROM"].astype(str) + ":" + df["POS"].astype(str) + ":" + df["REF"].astype(str) + ":" + df["ALT"].astype(str)

    # Genotype QC fields
    if has_gt:
        df["GT"] = df[gt_col].astype(str)
        df["GT_CLASS"] = df["GT"].map(parse_gt)
        df["DP"] = df[dp_col].map(as_int) if dp_col in df.columns else np.nan
        df["GQ"] = df[gq_col].map(as_int) if gq_col in df.columns else np.nan

        ref_alt_vaf = df[ad_col].map(parse_ad) if ad_col in df.columns else [(np.nan, np.nan, np.nan)] * len(df)
        ref_alt_vaf = np.array(list(ref_alt_vaf), dtype=float)
        df["AD_REF"] = ref_alt_vaf[:, 0]
        df["AD_ALT"] = ref_alt_vaf[:, 1]
        df["VAF"] = ref_alt_vaf[:, 2]
    else:
        df["GT"] = "."
        df["GT_CLASS"] = "MISSING"
        df["DP"] = np.nan
        df["GQ"] = np.nan
        df["AD_REF"] = np.nan
        df["AD_ALT"] = np.nan
        df["VAF"] = np.nan

    # Frequency filter
    freq_pass = df["AF_max"].isna() | (df["AF_max"] <= args.maf)

    # Genotype filter: keep only non-reference calls by default
    if has_gt:
        gt_nonref = df["GT_CLASS"].isin(["HET", "HOM_ALT"])
    else:
        gt_nonref = pd.Series(True, index=df.index) if args.allow_missing_genotype else pd.Series(False, index=df.index)

    # QC: DP, GQ, ALT depth, VAF
    qc_pass = pd.Series(True, index=df.index)

    if has_gt:
        qc_pass &= df["DP"].isna() | (df["DP"] >= args.min_dp)
        qc_pass &= df["GQ"].isna() | (df["GQ"] >= args.min_gq)
        qc_pass &= df["AD_ALT"].isna() | (df["AD_ALT"] >= args.min_alt_depth)

        # VAF constraints by genotype class
        het_mask = df["GT_CLASS"].eq("HET")
        hom_mask = df["GT_CLASS"].eq("HOM_ALT")

        # only apply when VAF is available
        qc_pass &= (~het_mask) | df["VAF"].isna() | ((df["VAF"] >= args.het_vaf_min) & (df["VAF"] <= args.het_vaf_max))
        qc_pass &= (~hom_mask) | df["VAF"].isna() | (df["VAF"] >= args.hom_vaf_min)

    # Functional tiering
    tiers = []
    reasons = []
    for _, row in df.iterrows():
        t, r = tier_and_reason(row, args)
        tiers.append(t)
        reasons.append(";".join(r))
    df["TIER"] = tiers
    df["REASON"] = reasons

    # Base keep mask
    base_keep = freq_pass & gt_nonref & qc_pass & df["TIER"].isin(["A", "B", "C"])

    kept = df[base_keep].copy()
    filtered = df[~base_keep].copy()

    # Inheritance model enforcement
    if args.mode == "AR" and has_gt:
        # Keep HOM_ALT directly
        hom_keep = kept["GT_CLASS"].eq("HOM_ALT")

        # Compound-het: genes with >=2 HET variants
        het = kept[kept["GT_CLASS"].eq("HET")].copy()
        comphet_genes = set(het[args.col_gene].value_counts()[lambda s: s >= 2].index)

        comphet_keep = kept[args.col_gene].isin(comphet_genes) & kept["GT_CLASS"].eq("HET")

        ar_keep = hom_keep | comphet_keep
        filtered = pd.concat([filtered, kept[~ar_keep]], ignore_index=True)
        kept = kept[ar_keep].copy()

        if args.emit_comphet:
            comphet_summary = (
                kept[kept["GT_CLASS"].eq("HET")]
                .groupby(args.col_gene)
                .agg(n_het=("VARIANT_KEY", "count"),
                     variants=("VARIANT_KEY", lambda x: ";".join(x.tolist())))
                .reset_index()
                .sort_values(["n_het"], ascending=False)
            )
            comphet_summary.to_csv(f"{args.out_prefix}.comphet_genes.tsv", sep="\t", index=False)

    # Output columns (keep your key fields first)
    preferred_cols = [
        "TIER", "REASON", "VARIANT_KEY",
        "CHROM", "POS", "REF", "ALT", "ID",
        args.col_gene, "ANN[0].FEATUREID", args.col_effect, args.col_impact, args.col_biotype,
        args.col_hgvsc, args.col_hgvsp,
        "AF_max",
        "gnomAD_cadd_phred", "REVEL_max", "SpliceAI_DS_max", "am_pathogenicity", "am_class",
        "GT", "GT_CLASS", "DP", "GQ", "AD_REF", "AD_ALT", "VAF"
    ]
    # Append all AF columns present
    af_cols = [c for c in df.columns if c.startswith("gnomAD_AF_") or c.startswith("WBBC_AF")]
    preferred_cols = [c for c in preferred_cols if c in df.columns] + [c for c in af_cols if c in df.columns]

    kept.sort_values(["TIER", "AF_max"], ascending=[True, True], inplace=True)

    out_prefix = Path(args.out_prefix)
    kept.to_csv(f"{out_prefix}.candidates.tsv", sep="\t", index=False, columns=preferred_cols)
    filtered.to_csv(f"{out_prefix}.filtered_out.tsv", sep="\t", index=False, columns=preferred_cols)

    # Summary
    summary = []
    summary.append(f"Input rows: {len(df)}")
    summary.append(f"Kept rows: {len(kept)}")
    summary.append(f"Filtered rows: {len(filtered)}")
    summary.append(f"Mode: {args.mode}")
    summary.append(f"MAF cutoff (AF_max): {args.maf}")
    summary.append(f"QC: DP>={args.min_dp}, GQ>={args.min_gq}, ALT_DP>={args.min_alt_depth}, "
                   f"HET_VAF in [{args.het_vaf_min},{args.het_vaf_max}], HOM_VAF>={args.hom_vaf_min}")
    summary.append(f"Predictors: SpliceAI support/strong={args.spliceai_support}/{args.spliceai_strong}; "
                   f"REVEL support/strong={args.revel_support}/{args.revel_strong}; "
                   f"CADD support/strong={args.cadd_support}/{args.cadd_strong}; "
                   f"AlphaMissense LP>={args.alphamissense_lp}")
    Path(f"{out_prefix}.summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
