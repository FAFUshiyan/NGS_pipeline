
suppressPackageStartupMessages({
  library(tidyverse)
  library(qqman)
  library(biomaRt)
})

# ===== 配置 =====
FILE        <- "~/Downloads/skat_SKATO_robust.all_masks.tsv"   # ← 改成你的SKATO结果路径
MASK_VALUE  <- "misd"                               # 只画 missense-damaging
GENOME      <- "GRCh38"                             # 或 "GRCh37"
OUT_PDF     <- "manhattan_misd.pdf"
LABEL_GENES <- c("PLAAT4")                             # 想特别标注的基因，可增删

# ===== 读取并识别列 =====
dat <- read.table(FILE, header = TRUE, sep = ifelse(grepl("\\.csv$", FILE, TRUE), ",", "\t"),
                  check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)

pick_first <- function(x) { y <- intersect(names(dat), x); if (length(y)) y[1] else NA_character_ }

col_gene <- pick_first(c("SetID","gene","Gene","GENE","SYMBOL","HGNC","GeneSymbol"))
col_mask <- pick_first(c("mask","Mask","MASK","category"))
col_p    <- pick_first(c("p_skato_robust","p_skato","p_acato","P.value","p.value","P","p"))
col_chr  <- pick_first(c("CHR","chr","chrom","Chromosome","chromosome"))
col_bp   <- pick_first(c("BP","bp","POS","pos","Start","start","begin","start_position"))

stopifnot(!is.na(col_gene), !is.na(col_p))

res <- dat %>%
  { if (!is.na(col_mask)) dplyr::filter(., .data[[col_mask]] == MASK_VALUE) else . } %>%
  transmute(
    gene = .data[[col_gene]],
    P    = suppressWarnings(as.numeric(.data[[col_p]])),
    CHR  = if (!is.na(col_chr)) .data[[col_chr]] else NA,
    BP   = if (!is.na(col_bp))  .data[[col_bp]]  else NA
  ) %>%
  distinct(gene, .keep_all = TRUE)

# ===== 若无CHR/BP，用 biomaRt 取基因坐标（取中点） =====
need_anno <- any(is.na(res$CHR)) || any(is.na(res$BP))
if (need_anno) {
  message("没有提供CHR/BP，使用 biomaRt 注释基因坐标…")
  if (GENOME == "GRCh37") {
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                       host = "grch37.ensembl.org")
  } else {
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  }
  ids <- res$gene
  id_type <- if (all(grepl("^ENSG", ids, ignore.case = TRUE))) "ensembl_gene_id" else "hgnc_symbol"
  
  ann <- getBM(
    attributes = c(id_type, "chromosome_name","start_position","end_position"),
    filters    = id_type,
    values     = unique(ids),
    mart       = mart
  ) %>%
    mutate(BP = floor((start_position + end_position)/2)) %>%
    transmute(gene_id = .data[[id_type]], CHR = chromosome_name, BP)
  
  res <- res %>%
    left_join(ann, by = c("gene" = "gene_id")) %>%
    mutate(
      CHR = coalesce(CHR.x, CHR.y),
      BP  = coalesce(BP.x,  BP.y)
    ) %>%
    select(gene, P, CHR, BP)
}

# ===== 清洗并转成 qqman 需要的格式 =====
res <- res %>%
  mutate(
    CHR = toupper(as.character(CHR)),
    CHR = recode(CHR, "X"="23", "Y"="24", "MT"="25", "M"="25"),
    CHR = suppressWarnings(as.integer(CHR)),
    BP  = suppressWarnings(as.integer(BP))
  ) %>%
  filter(!is.na(CHR), !is.na(BP), !is.na(P), P > 0 & P <= 1) %>%
  rename(SNP = gene)

# ===== 计算按“基因数”Bonferroni 阈值并作图 =====
bonf <- 0.05 / dplyr::n_distinct(res$SNP)
ylim_top <- max(-log10(res$P), na.rm = TRUE) + 1
###下面这行代码与是否要X，Y染色体相关
res <- res %>% dplyr::filter(CHR %in% 1:22)
pdf(OUT_PDF, width = 12, height = 5)
manhattan(
  res, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
  genomewideline = FALSE, suggestiveline = FALSE,
  cex = 0.6, cex.axis = 0.9, ylim = c(0, ylim_top),
  highlight = intersect(LABEL_GENES, res$SNP)
)
abline(h = -log10(bonf), lty = 2, col = "red")  # 基因级Bonferroni阈值
mtext(sprintf("Bonferroni (0.05 / %d genes) = %.2g", dplyr::n_distinct(res$SNP), bonf),
      side = 3, adj = 1, cex = 0.8)
dev.off()

cat("Done! 输出：", OUT_PDF, "\n")
