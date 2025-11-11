# ========= 依赖 =========
# install.packages(c("SKAT","data.table"))  # 如未安装先装
library(SKAT)
library(data.table)

# ========= 路径参数（按需改） =========
bed_prefix <- "test.bed1u"            # 你的 PLINK 前缀：.bed/.bim/.fam
pheno_file <- "pheno_covars.tsv"      # 你刚上传并确认好的文件



setids <- c(
  lof    = "setid_lof.txt",
  misd   = "setid_misd.txt",
  splice = "setid_splice.txt",
  merged = "setid_merged.txt"
)

# SSD 输出目录
ssd_dir <- "ssd_out"
dir.create(ssd_dir, showWarnings = FALSE)
gen_ssd_if_needed <- function(setid_file, ssd_prefix) {
  file_ssd  <- file.path(ssd_dir, paste0(ssd_prefix, ".SSD"))
  file_info <- file.path(ssd_dir, paste0(ssd_prefix, ".SSD.info"))
  if (!file.exists(file_ssd) || !file.exists(file_info)) {
    message("Generating SSD for: ", setid_file)
    Generate_SSD_SetID(
      File.Bed   = bed_path,
      File.Bim   = bim_path,
      File.Fam   = fam_path,
      File.SetID = setid_file,
      File.SSD   = file_ssd,
      File.Info  = file_info,
      Is.FlipGenotype = TRUE
    )
  }
  return(list(ssd=file_ssd, info=file_info))
}

ssd_list <- lapply(names(setids), function(mask) {
  gen_ssd_if_needed(setids[[mask]], paste0("mask_", mask))
})
names(ssd_list) <- names(setids)
# 如果你已经生成了 SSD，直接写这里；否则在下面的“可选：自动生成SSD”里按需生成
ssd <- list(
  lof    = list(SSD="ssd_out/mask_lof.SSD",    INFO="ssd_out/mask_lof.SSD.info"),
  misd   = list(SSD="ssd_out/mask_misd.SSD",   INFO="ssd_out/mask_misd.SSD.info"),
  splice = list(SSD="ssd_out/mask_splice.SSD", INFO="ssd_out/mask_splice.SSD.info"),
  merged = list(SSD="ssd_out/mask_merged.SSD", INFO="ssd_out/mask_merged.SSD.info")
)

# ========= 1) 读取 pheno+协变量，并与 .fam 对齐 =========
fam_path <- paste0(bed_prefix, ".fam")
bim_path <- paste0(bed_prefix, ".bim")
bed_path <- paste0(bed_prefix, ".bed")

ph <- fread(pheno_file, sep="\t", header=TRUE, data.table=FALSE)
names(ph) <- trimws(sub("^#","", names(ph)))
stopifnot(all(c("IID","pheno","sex") %in% names(ph)))

# 自动识别 PCs（存在多少用多少）
pc_cols <- grep("^PC\\d+$", names(ph), value=TRUE)
stopifnot(length(pc_cols) >= 1)

# 类型统一
ph$IID   <- as.character(ph$IID)
ph$pheno <- as.numeric(ph$pheno)  # 0/1
ph$sex   <- as.numeric(ph$sex)    # 1/2（也可 factor）

# fam 命名
fam <- fread(fam_path, header=FALSE, data.table=FALSE)
colnames(fam)[1:6] <- c("FID","IID","PID","MID","SEX_FAM","PHENO_FAM")
fam$IID <- as.character(fam$IID)

# 按 fam 的 IID 顺序左连接（与基因型顺序一致）
setDT(ph); setDT(fam)
ph_aln <- ph[fam[,.(IID)], on=.(IID)]
stopifnot(!any(is.na(ph_aln$IID)))
stopifnot(anyDuplicated(ph_aln$IID) == 0)

# 去除自变量/因变量缺失
model_vars <- c("pheno","sex", pc_cols)
ok <- stats::complete.cases(ph_aln[, model_vars, with=FALSE])
ph_model <- as.data.frame(ph_aln[ok])

cat(sprintf("[INFO] 样本量用于空模型：%d\n", nrow(ph_model)))

# ========= 2) 建立空模型（二分类；小样本会自动小样本校正） =========
form <- as.formula(paste("pheno ~ sex +", paste(pc_cols, collapse=" + ")))
obj <- SKAT_Null_Model(form, out_type="D", data=ph_model)
# 文档：先用 SKAT_Null_Model 得到空模型，再把对象传给 SKAT/…SSD.All 家族。:contentReference[oaicite:1]{index=1}

# ========= （可选）如果还没生成 SSD：按 SetID 现生成 =========
# 仅在你还没 .SSD/.info 时启用；每个 setid_*.txt 两列：Gene <tab> VariantID
# 且第二列“VariantID”必须与 .bim 第二列完全一致（如 chr:pos:REF:ALT）。
#Generate_SSD_SetID 用法见文档。:contentReference[oaicite:2]{index=2}
#dir.create("ssd_out", showWarnings=FALSE)
#Generate_SSD_SetID(bed_path, bim_path, fam_path,
#                   "setid_lof.txt",    "ssd_out/mask_lof.SSD",    "ssd_out/mask_lof.SSD.info",    Is.FlipGenotype=TRUE)
#Generate_SSD_SetID(bed_path, bim_path, fam_path,
#                   "setid_misd.txt",   "ssd_out/mask_misd.SSD",   "ssd_out/mask_misd.SSD.info",   Is.FlipGenotype=TRUE)
#Generate_SSD_SetID(bed_path, bim_path, fam_path,
#                   "setid_splice.txt", "ssd_out/mask_splice.SSD", "ssd_out/mask_splice.SSD.info", Is.FlipGenotype=TRUE)
#Generate_SSD_SetID(bed_path, bim_path, fam_path,
                    "setid_merged.txt", "ssd_out/mask_merged.SSD", "ssd_out/mask_merged.SSD.info", Is.FlipGenotype=TRUE)

# ========= 3) 跑稳健 SKAT-O（逐面具） =========
run_mask <- function(mask_name, s) {
  if (!file.exists(s$SSD) || !file.exists(s$INFO)) {
    warning("找不到 SSD：", mask_name); return(NULL)
  }
  SSD.INFO <- Open_SSD(s$SSD, s$INFO)        # 打开 SSD（用完要 Close_SSD）:contentReference[oaicite:3]{index=3}
  on.exit(Close_SSD(), add=TRUE)

  # 稳健二分类接口，method="SKATO" → SKAT-O；权重 Beta(1,25) 是常用默认。:contentReference[oaicite:4]{index=4}
  out <- SKATBinary_Robust.SSD.All(
    SSD.INFO, obj,
    method="SKATO",
    kernel="linear.weighted",
    weights.beta=c(1,25)
  )
  res <- out$results
  if (!is.null(res) && nrow(res) > 0) {
    res$FDR_BH <- p.adjust(res$P.value, method="BH")
    res$mask   <- mask_name
  }
  res
}

res_list <- lapply(names(ssd), function(m) run_mask(m, ssd[[m]]))
res <- do.call(rbind, Filter(Negate(is.null), res_list))
fwrite(res, file="skat_SKATO_robust.all_masks.tsv", sep="\t")
cat("[DONE] 结果已写入 skat_SKATO_robust.all_masks.tsv\n")
