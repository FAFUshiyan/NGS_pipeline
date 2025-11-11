# =========================================================
# ACAT-O over SKAT-O results — FULL DEBUG VERSION
# - 仅保留: lof / misd / splice；显式排除 merged
# - 长表或宽表皆可；每一步输出到 acato_debug/
# - ACAT-O: sumFREGAT::ACATO(p) ；等权ACAT: acat_cct() 纯R实现
# =========================================================
suppressPackageStartupMessages({
  library(data.table)
})

IN  <- "~/Downloads/skat_SKATO_robust.all_masks.tsv"  # ← 改成你的输入
OUT <- "genes_ACATO_from_SKATO.tsv"
DBG <- "acato_debug"
dir.create(DBG, showWarnings = FALSE, recursive = TRUE)

# ---------- helpers ----------
writetsv <- function(x, path, nmax_head = 50L){
  fp <- file.path(DBG, path)
  fwrite(x, fp, sep = "\t")
  hh <- tryCatch(head(x, nmax_head), error=function(e) NULL)
  if (!is.null(hh)) fwrite(hh, sub("\\.tsv$", ".head.tsv", fp), sep = "\t")
}
writetxt <- function(txt, path){
  cat(paste0(txt, collapse = "\n"), file = file.path(DBG, path))
}
logmsg <- function(...) cat(sprintf("[acato] %s\n", sprintf(...)))

options(error = NULL)  # 防止 options(error=recover) 进入 Browse[...]

# 0) 读入
DT <- fread(IN, na.strings = c("NA","NaN",""))
writetsv(DT, "00_raw.tsv"); logmsg("read %d rows, %d cols", nrow(DT), ncol(DT))

# 1) 列名识别
cn  <- names(DT); low <- tolower(cn)
gene_col <- cn[match(TRUE, low %in% c("gene","setid"))]
mask_col <- cn[match(TRUE, low %in% c("mask","masks","class"))]
p_candidates <- c("p","p.value","p_value","pvalue","p.skato","p_skato","skato_p","pval","pval_skato")
p_col <- cn[match(TRUE, low %in% p_candidates)]

writetxt(c(
  sprintf("gene_col_detected: %s", ifelse(is.na(gene_col), "NA", gene_col)),
  sprintf("mask_col_detected: %s", ifelse(is.na(mask_col), "NA", mask_col)),
  sprintf("p_col_detected   : %s", ifelse(is.na(p_col), "NA", p_col))
), "01_detected_columns.txt")

if (is.na(gene_col)) stop("未找到基因列（Gene/SetID）")
DT[, (gene_col) := trimws(as.character(get(gene_col)))]

# 允许&保留的面具（严控只用三类）
allow_masks <- c("lof","misd","splice")

# 2) 统一为【宽表】，并显式排除 merged
if (!is.na(mask_col) && !is.na(p_col)) {
  # ---- 长表路径 ----
  DT[, (mask_col) := tolower(trimws(as.character(get(mask_col))))]
  DT[, (p_col)    := suppressWarnings(as.numeric(get(p_col)))]
  DT <- DT[is.finite(get(p_col)) & get(p_col) > 0 & get(p_col) < 1]
  writetsv(DT, "02_long_clean.tsv")

  # 2a) 排除 merged；仅保留 allow_masks
  DT_excl_merged <- DT[tolower(get(mask_col)) == "merged"]
  writetsv(DT_excl_merged, "02a_long_excluded_merged.tsv")
  DT <- DT[tolower(get(mask_col)) != "merged"]
  DT <- DT[tolower(get(mask_col)) %in% allow_masks]
  writetsv(DT, "02b_long_after_excluding_merged_and_filter_masks.tsv")

  # (gene,mask) 多行 → 取最小 p（如需先对同一(基因,面具)内部ACAT再合，也可扩展）
  DTmin <- DT[, .(p = min(get(p_col), na.rm = TRUE)), by = c(gene_col, mask_col)]
  setnames(DTmin, c(gene_col, mask_col), c("Gene",".mask"))
  writetsv(DTmin, "03_long_collapsed_min.tsv")

  # 宽表
  wide <- dcast(DTmin, Gene ~ .mask, value.var = "p")
  writetsv(wide, "04_wide_raw.tsv")
  # 统一列名 p_<mask>
  mcols <- setdiff(names(wide), "Gene")
  setnames(wide, mcols, paste0("p_", mcols), skip_absent = TRUE)
  gene_col <- "Gene"

} else {
  # ---- 宽表路径 ----
  wide <- copy(DT)
  setnames(wide, gene_col, "Gene")
  gene_col <- "Gene"
  writetsv(wide, "02B_wide_input.tsv")

  # 删除 p_merged* 列
  merged_cols <- grep("^p[_.]?merged", tolower(names(wide)), value = TRUE)
  if (length(merged_cols)) {
    writetxt(paste(merged_cols, collapse = "\n"), "02B_removed_merged_columns.txt")
    wide[, (merged_cols) := NULL]
  } else {
    writetxt("no p_merged-like columns detected", "02B_removed_merged_columns.txt")
  }
  # 只保留 allow_masks
  keep_regex <- sprintf("^p[_.]?(%s)$", paste(allow_masks, collapse = "|"))
  keep_cols  <- grep(keep_regex, tolower(names(wide)), value = TRUE)
  drop_cols  <- grep("^p[_.]?", tolower(names(wide)), value = TRUE)
  drop_cols  <- setdiff(drop_cols, keep_cols)
  if (length(drop_cols)) {
    writetxt(paste(drop_cols, collapse = "\n"), "02B_dropped_non_allowed_masks.txt")
    wide[, (drop_cols) := NULL]
  }
  writetsv(wide, "02C_wide_after_mask_filter.tsv")
}

# 3) 面具列数值化（严格只用 lof/misd/splice）
mask_regex <- "^p[_.]?(lof|misd|splice)$"
found_cols <- grep(mask_regex, tolower(names(wide)), value = TRUE)
std_cols   <- sub("^p[_.]?", "p_", tolower(found_cols))
setnames(wide, found_cols, std_cols, skip_absent = TRUE)

mask_cols <- intersect(names(wide), std_cols)
if (length(mask_cols) == 0L) stop("未识别到用于合并的面具列（仅允许 p_lof/p_misd/p_splice）")

for (c in mask_cols) wide[, (c) := suppressWarnings(as.numeric(get(c)))]
writetsv(wide[, c(gene_col, mask_cols), with = FALSE], "05_wide_std_numeric.tsv")
writetxt(paste(mask_cols, collapse = "\n"), "05_mask_cols_used.txt")

# 4) 删除“所有面具均 NA”的基因（.SD 仅在 j 中使用）
keep_vec <- wide[, rowSums(!is.na(.SD)) > 0, .SDcols = mask_cols]
wide_keepinfo <- data.table(Gene = wide[[gene_col]],
                            n_nonNA = wide[, rowSums(!is.na(.SD)), .SDcols = mask_cols],
                            keep = keep_vec)
writetsv(wide_keepinfo, "06_keep_vector.tsv")
wide <- wide[keep_vec]
writetsv(wide, "07_wide_kept.tsv")

# 5) ACAT-O 与等权 ACAT
# 5.1 ACAT-O: 官方实现（输入向量p→输出组合p）
if (!requireNamespace("sumFREGAT", quietly = TRUE)) {
  stop("需要 sumFREGAT 包：install.packages('sumFREGAT')")
}
ACATO_safe <- function(p){
  p <- as.numeric(p); p <- p[is.finite(p) & !is.na(p)]
  if (!length(p)) return(NA_real_)
  if (length(p) == 1L) return(max(min(p, 1), .Machine$double.xmin))
  sumFREGAT::ACATO(p)  # ACAT-O：向量→单个组合p
}

# 5.2 等权 ACAT：纯R实现（柯西变换→均值→反变换）
acat_cct <- function(p){
  p <- as.numeric(p); p <- p[is.finite(p) & !is.na(p)]
  if (!length(p)) return(NA_real_)
  p[p < 1e-15]     <- 1e-15
  p[p > 1 - 1e-15] <- 1 - 1e-15
  Tsum <- mean(tan((0.5 - p) * pi))
  pc   <- 0.5 - atan(Tsum) / pi
  max(min(pc, 1.0), .Machine$double.xmin)
}

# 为便于核查，输出“长格式”的入参
long_input <- melt(wide, id.vars = gene_col, measure.vars = mask_cols,
                   variable.name = "mask", value.name = "p")
long_input[, mask := sub("^p_", "", mask)]
writetsv(long_input, "08_acat_inputs_long.tsv")

# 逐行计算
wide[, p_acato := apply(.SD, 1, ACATO_safe), .SDcols = mask_cols]
wide[, p_acat  := apply(.SD, 1, acat_cct  ), .SDcols = mask_cols]
writetsv(wide[, c(gene_col, mask_cols, "p_acato", "p_acat"), with = FALSE], "09_with_acat_cols.tsv")

# 6) 多重校正 + 排序 + 导出
wide[, FDR_BH := p.adjust(p_acato, method = "BH")]
wide[, n_masks := rowSums(!is.na(.SD)), .SDcols = mask_cols]
data.table::setorderv(wide, "p_acato", na.last = TRUE)
writetsv(wide, "10_final_ranked.tsv")

fwrite(wide, OUT, sep = "\t")
cat(sprintf("Done. 写出主结果：%s\n中间文件目录：%s\n", OUT, normalizePath(DBG, winslash = '/')))
