#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build SKAT SetID files from a TSV with SnpEff/SnpSift columns,
matching PLINK .bim SNP ID style.

- ID styles:
    colon:  chrom:pos:REF:ALT
    comma:  chrom:posREF,ALT     <-- your .bim style (e.g., 1:13273G,C)
- chr prefix:
    yes/no/auto

Outputs: setid_lof.txt / setid_misd.txt / setid_splice.txt / setid_merged.txt
"""

import argparse, sys, re
from collections import defaultdict

ALIASES = {
    "chrom": ["CHROM", "#CHROM", "Chrom", "chrom"],
    "pos":   ["POS", "Position", "pos"],
    "ref":   ["REF", "ref"],
    "alt":   ["ALT", "alt"],
    "effect":["ANN[*].EFFECT", "EFFECT", "Annotation", "Consequence"],
    "gene":  ["ANN[*].GENE", "GENE", "Gene_Name", "Gene", "SYMBOL"],
    "biotype":["ANN[*].BIOTYPE", "BIOTYPE", "Transcript_biotype", "Biotype"],
    "ann":   ["ANN", "ANN_FIRST", "ANN1"]
}

def find_col(cands, header):
    idx = {c:i for i,c in enumerate(header)}
    for c in cands:
        if c in idx: return c, idx[c]
    return None, None

def load_bim_ids(bim):
    ids=set()
    if not bim: return ids
    with open(bim) as f:
        for line in f:
            p=line.rstrip("\n").split()
            if len(p)>=2: ids.add(p[1])
    return ids

def make_id(chrom,pos,ref,alt,style,with_chr):
    # 统一/去除 chr 前缀
    base = ("chr"+chrom.lstrip("chr")) if with_chr else chrom.lstrip("chr")
    if style=="comma":
        return f"{base}:{pos}{ref},{alt}"
    else:
        return f"{base}:{pos}:{ref}:{alt}"

def float0(x):
    try:
        if x in (None,"",".","NA"): return 0.0
        return float(x)
    except: return 0.0

def parse_first_ann(ann_raw):
    first = ann_raw.split(",")[0]
    parts = first.split("|")
    eff   = parts[1] if len(parts)>1 else "."
    imp   = parts[2] if len(parts)>2 else "."
    gene  = parts[3] if len(parts)>3 else "."
    biot  = parts[7] if len(parts)>7 else "."
    return eff, imp, gene, biot

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tsv", help="注释 TSV（含 CHROM/POS/REF/ALT 和效应/打分/频率列）")
    ap.add_argument("--bim", help="cohort.bim（用于ID对齐与门禁）")
    ap.add_argument("--maf", type=float, default=0.01, help="maxAF 阈值（默认 0.01）")
    ap.add_argument("--revel", type=float, default=0.5)
    ap.add_argument("--alphamiss", type=float, default=0.564)
    ap.add_argument("--cadd", type=float, default=20.0)
    ap.add_argument("--spliceai", type=float, default=0.5)
    ap.add_argument("--id-style", choices=["auto","comma","colon"], default="auto",
                    help="变异ID风格；你的 .bim 是 comma")
    ap.add_argument("--force-chr", choices=["auto","yes","no"], default="auto",
                    help="是否带 chr 前缀；你的 .bim 用 'no'")
    args = ap.parse_args()

    # 读表头
    with open(args.tsv, encoding="utf-8") as f:
        header=f.readline().rstrip("\n").split("\t")

    # 关键列
    cname, ci = find_col(ALIASES["chrom"], header)
    pname, pi = find_col(ALIASES["pos"], header)
    rname, ri = find_col(ALIASES["ref"], header)
    aname, ai = find_col(ALIASES["alt"], header)
    if None in (ci,pi,ri,ai):
        sys.exit("[ERROR] 缺少 CHROM/POS/REF/ALT 列")

    ename, ei = find_col(ALIASES["effect"], header)
    gname, gi = find_col(ALIASES["gene"], header)
    bname, bi = find_col(ALIASES["biotype"], header)
    ann_name, ann_i = find_col(ALIASES["ann"], header)

    have_split = (ei is not None and gi is not None and bname is not None and bi is not None)
    have_ann = (ann_i is not None)
    if not have_split and not have_ann:
        sys.exit("[ERROR] 找不到 EFFECT/GENE/BIOTYPE 拆分列，也没有 ANN 列")

    # 分数/频率列
    col_idx = {c:i for i,c in enumerate(header)}
    af_cols = [c for c in header if c.startswith("gnomAD_AF_")] + (["WBBC_AF"] if "WBBC_AF" in header else [])
    def GET(row, key): return row[col_idx[key]] if key in col_idx and col_idx[key] < len(row) else "."

    # .bim ID 集合
    bim_ids = load_bim_ids(args.bim)

    # 决定 chr 前缀和风格：你的 .bim 是 “comma + no chr”
    if args.id_style == "auto":
        id_style = "comma" if any("," in x and ":" in x and not x.startswith("chr") for x in list(bim_ids)[:1000]) else "colon"
    else:
        id_style = args.id_style
    if args.force_chr == "auto":
        with_chr = False if any(not x.startswith("chr") for x in list(bim_ids)[:1000]) else True
    else:
        with_chr = (args.force_chr == "yes")

    # 容器
    buckets = {"lof":defaultdict(set), "misd":defaultdict(set), "splice":defaultdict(set)}

    # 遍历
    seen = 0; hits = 0
    with open(args.tsv, encoding="utf-8") as f:
        next(f)
        for line in f:
            if not line.strip(): continue
            row=line.rstrip("\n").split("\t")
            chrom,pos,ref,alt = row[ci], row[pi], row[ri], row[ai]
            vid = make_id(chrom,pos,ref,alt,id_style,with_chr)

            seen += 1
            if bim_ids:
                if vid not in bim_ids: 
                    continue
                hits += 1

            if have_split:
                effect = row[ei]
                gene   = row[gi]
                biotype= row[bi]
            else:
                effect, impact, gene, biotype = parse_first_ann(row[ann_i])

            if not gene or gene==".":
                continue

            # maxAF
            af_max = max([float0(GET(row,c)) for c in af_cols], default=0.0)
            if af_max >= args.maf: 
                continue

            # scores
            revel = float0(GET(row,"REVEL_max"))
            cadd  = float0(GET(row,"gnomAD_cadd_phred")) if "gnomAD_cadd_phred" in header else float0(GET(row,"CADD_PHRED"))
            ams   = float0(GET(row,"am_pathogenicity"))
            sai   = float0(GET(row,"SpliceAI_DS_max"))

            is_pc = ("protein_coding" in (biotype or ""))
            is_lof = bool(re.search(r"(splice_acceptor_variant|splice_donor_variant|stop_gained|frameshift_variant|start_lost|stop_lost|transcript_ablation)", effect))
            is_mis = (effect == "missense_variant")
            is_sp  = bool(re.search(r"(splice_acceptor_variant|splice_donor_variant)", effect)) or (sai >= args.spliceai)

            if is_pc and is_lof: buckets["lof"][gene].add(vid)
            if is_mis and ((revel>=args.revel) or (ams>=args.alphamiss) or (cadd>=args.cadd)): buckets["misd"][gene].add(vid)
            if is_sp: buckets["splice"][gene].add(vid)

    # 输出
    def dump(path, d):
        with open(path,"w") as out:
            for g in sorted(d.keys()):
                for v in sorted(d[g]):
                    out.write(f"{g}\t{v}\n")

    dump("setid_lof.txt",    buckets["lof"])
    dump("setid_misd.txt",   buckets["misd"])
    dump("setid_splice.txt", buckets["splice"])

    merged = defaultdict(set)
    for k in buckets:
        for g, vs in buckets[k].items():
            merged[g].update(vs)
    dump("setid_merged.txt", merged)

    # 汇总/命中率
    def sz(d): 
        return len(d), sum(len(vs) for vs in d.values())
    print(f"[ID style] {id_style}  [chr prefix] {'yes' if with_chr else 'no'}")
    if bim_ids:
        rate = (hits/seen*100) if seen else 0
        print(f"[.bim gate] seen={seen}  matched={hits}  match_rate={rate:.1f}%")
    for name, d in [("LoF",buckets["lof"]),("MisD",buckets["misd"]),("Splice",buckets["splice"]),("Merged",merged)]:
        gcnt, vcnt = sz(d)
        print(f"[{name}] genes={gcnt} variants={vcnt}")

if __name__ == "__main__":
    main()

