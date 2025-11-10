# 必要包
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)  # 用于标签自动避让
})

# ===== 1) 读入并整理 =====
# 将路径改成你的文件位置；若与脚本同目录，直接用文件名即可
infile <- "~/Downloads/lipid.KS.tsv"

dat_raw <- read_tsv(infile, show_col_types = FALSE)

# 适配你附件中的列名：Component Name / log2FC / (-log10 Pvalue
dat <- dat_raw %>%
  rename(
    name   = `Component Name`,
    log2FC = `log2FC`,
    mlog10p = `-log10 Pvalue`
  ) %>%
  mutate(
    sig = mlog10p > 1 & abs(log2FC) > 0.5,
    direction = case_when(
      sig & log2FC >  0.5 ~ "Up",
      sig & log2FC < -0.5 ~ "Down",
      TRUE                ~ "NS"
    )
  )

# 显著性计数（可查看或写出）
count_tab <- dat %>% count(direction)
print(count_tab)

# ===== 2) 可选：挑选需要标注的点（避免过度拥挤）=====
# 这里选择显著点里 -log10(p) 排名前15 的做标签；如不想标注可将 p_label 设为 NULL
p_label <- dat %>%
  filter(sig) %>%
  arrange(desc(mlog10p)) %>%
  slice_head(n = 15)

# ===== 3) 绘图 =====
p <- ggplot(dat, aes(x = log2FC, y = mlog10p, color = direction)) +
  geom_point(alpha = 0.85, size = 2) +
  # 阈值线
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5) +
  # 颜色可按需调整
  scale_color_manual(values = c(Down = "#1f77b4", NS = "grey70", Up = "#d62728")) +
  labs(
    title = "Volcano plot",
    x = "log2 fold change",
    y = "-log10(p-value)",
    color = "Status"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )+coord_cartesian(xlim = c(-3, 3)) 

# 可选标签（显著且Top 15）
if (!is.null(p_label) && nrow(p_label) > 0) {
  p <- p + ggrepel::geom_text_repel(
    data = p_label,
    aes(label = name),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0,
    box.padding = 0.3
  )
}
p+coord_cartesian(xlim = c(-3, 3)) 
print(p)

# ===== 4) 导出 =====
ggsave("volcano.png", p, width = 6, height = 5, dpi = 300)
ggsave("volcano.pdf", p, width = 6, height = 5)
