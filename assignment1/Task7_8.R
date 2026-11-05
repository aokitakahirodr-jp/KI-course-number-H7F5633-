## 7.1 magic_guys.csv

```{r}
magic <- read.csv("~/KI-course-number-H7F5633-/magic_guys.csv")
str(magic)
table(magic$species)
```

### a. Histograms

```{r, fig.height = 7}
par(mfrow = c(2, 2))                       # 2x2 panel: three break settings + legend
for (b in c(5, 10, 20)) {                  # overlay jedi (blue) and sith (red) per break setting
  hist(magic$length[magic$species == "jedi"], breaks = b, col = rgb(0, 0, 1, 0.4),
       xlim = range(magic$length), main = paste("hist(), breaks =", b), xlab = "Body height")
  hist(magic$length[magic$species == "sith"], breaks = b, col = rgb(1, 0, 0, 0.4), add = TRUE)
}
plot.new(); legend("center", c("jedi", "sith"), fill = c(rgb(0, 0, 1, 0.4), rgb(1, 0, 0, 0.4)), bty = "n")
par(mfrow = c(1, 1))
```

```{r}
# overlapping histograms, binwidth 10
p_hist <- ggplot(magic, aes(length, fill = species)) +
  geom_histogram(binwidth = 10, alpha = 0.6, position = "identity", colour = "white") +
  scale_fill_manual(values = c(jedi = "steelblue", sith = "firebrick")) +
  labs(x = "Body height", y = "Count", title = "Body height by species") +
  theme_minimal()
p_hist

# alternative: facets and a different binwidth
ggplot(magic, aes(length, fill = species)) +
  geom_histogram(binwidth = 20, colour = "white") +
  facet_wrap(~ species, ncol = 1) +
  scale_fill_manual(values = c(jedi = "steelblue", sith = "firebrick")) +
  theme_bw() + theme(legend.position = "none")
```

Jedi are on average taller and more variable in height than Sith; the distributions overlap. Fewer breaks smooth the shape, more breaks reveal gaps caused by the small sample size (n = 50 per species).

### b. Boxplots

```{r}
# base boxplot
boxplot(length ~ species, data = magic, col = c("steelblue", "firebrick"),
        ylab = "Body height", main = "Body height by species")

# ggplot boxplot with individual points overlaid
p_box <- ggplot(magic, aes(species, length, fill = species)) +
  geom_boxplot(width = 0.5, outlier.shape = 21) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c(jedi = "steelblue", sith = "firebrick")) +
  labs(y = "Body height") +
  theme_minimal() + theme(legend.position = "none")
p_box
```

Boxplots confirm the higher median and wider interquartile range of Jedi heights.

### c. Saving plots

```{r}
# save each plot in raster (png) and vector (pdf, svg) formats
ggsave("figures/height_hist.png", p_hist, width = 7, height = 4.5, dpi = 300)
ggsave("figures/height_hist.pdf", p_hist, width = 7, height = 4.5)
ggsave("figures/height_hist.svg", p_hist, width = 7, height = 4.5)   # requires svglite
ggsave("figures/height_box.png", p_box, width = 5, height = 4.5, dpi = 300)
ggsave("figures/height_box.pdf", p_box, width = 5, height = 4.5)
ggsave("figures/height_box.svg", p_box, width = 5, height = 4.5)
list.files("figures")
```

- **PNG** (raster): slides, web pages, quick sharing; fixed resolution, small file for complex plots.
- **PDF** (vector): journal submission and print; scales without quality loss, embeds fonts, multi-page possible.
- **SVG** (vector): web embedding and further editing in Illustrator/Inkscape (e.g. figure assembly), editable text and shapes.

## 7.2 microarray_data.tab

```{r}
ma <- read.delim("~/KI-course-number-H7F5633-/microarray_data.tab")
dim(ma)   # a. rows = samples, columns = genes
```

The matrix has `r nrow(ma)` rows (samples) and `r ncol(ma)` columns (genes g1–g1000).

```{r}
# b. missing values per gene (column)
na_per_gene <- colSums(is.na(ma))
summary(na_per_gene)

# histogram of missing-value counts across genes
tibble(gene = names(na_per_gene), n_missing = na_per_gene) %>%
  ggplot(aes(n_missing)) +
  geom_histogram(binwidth = 10, fill = "grey40", colour = "white") +
  labs(x = "Missing values per gene", y = "Number of genes",
       title = "Distribution of missing values across 1000 genes") +
  theme_minimal()
```

```{r}
# c. genes exceeding X% missing
na_frac <- na_per_gene / nrow(ma)
sapply(c(0.1, 0.2, 0.5), function(x) sum(na_frac > x)) %>% setNames(c(">10%", ">20%", ">50%"))
head(names(na_frac)[na_frac > 0.5])   # example genes with >50 % missing
```

```{r}
# d. mean imputation per gene
ma_imputed <- ma %>% mutate(across(everything(), ~ replace(.x, is.na(.x), mean(.x, na.rm = TRUE))))
sum(is.na(ma_imputed))                 # remaining NAs
sum(na_per_gene == nrow(ma))           # genes with no observed value at all -> mean undefined
ma_imputed <- ma_imputed %>% select(where(~ !all(is.na(.x))))   # remove the all-NA genes
dim(ma_imputed); sum(is.na(ma_imputed))
```

`r sum(na_per_gene == nrow(ma))` genes have no observed value in any sample, so no mean exists and they are removed. Mean imputation preserves the per-gene mean but shrinks variance and ignores correlation structure; genes with a large fraction of missing values should be filtered rather than imputed.

## 7.3 CO2 dataset visualisation

```{r}
# one uptake curve per plant, coloured by treatment, panelled by origin
ggplot(CO2, aes(conc, uptake, colour = Treatment, group = Plant)) +
  geom_line(alpha = 0.6) + geom_point(size = 1.5) +
  facet_wrap(~ Type) +
  scale_colour_manual(values = c(nonchilled = "darkorange", chilled = "steelblue")) +
  labs(x = expression("Ambient CO"[2]*" (µL/L)"), y = expression("CO"[2]*" uptake (µmol/m"^2*"/s)"),
       title = "CO2 uptake curves per plant") +
  theme_bw()
```

Uptake saturates with increasing CO2 concentration in all plants. Quebec plants reach higher uptake than Mississippi plants, and chilling reduces uptake in both origins, with a much stronger effect in Mississippi plants.

# Task 8 — Tidybiology

```{r}
# devtools::install_github("hirscheylab/tidybiology")
suppressPackageStartupMessages(library(tidybiology))
data(chromosome); data(proteins)
```

```{r}
# a. mean/median/max per variable, reshaped to one row per variable
chromosome %>%
  summarise(across(c(variations, protein_codinggenes, mi_rna),
                   list(mean = mean, median = median, max = max), .names = "{.col}_{.fn}")) %>%
  pivot_longer(everything(), names_to = c("variable", "statistic"), names_pattern = "(.*)_(mean|median|max)") %>%
  pivot_wider(names_from = statistic, values_from = value)
```

```{r}
# b. chromosome length, ordered bar chart
ggplot(chromosome, aes(x = fct_reorder(id, basepairs, .desc = TRUE), y = basepairs / 1e6)) +
  geom_col(fill = "steelblue") +
  labs(x = "Chromosome", y = "Length (Mb)", title = "Human chromosome sizes") +
  theme_minimal()

# histogram of chromosome lengths
ggplot(chromosome, aes(basepairs / 1e6)) +
  geom_histogram(bins = 8, fill = "grey40", colour = "white") +
  labs(x = "Chromosome length (Mb)", y = "Count") + theme_minimal()
```

Sizes range from ~47 Mb (chr21) to ~249 Mb (chr1), roughly uniformly spread with chr1–chr12 and X above 130 Mb.

```{r, fig.height = 4}
# c. scatter + linear fit, points labelled by chromosome
p_pc <- ggplot(chromosome, aes(basepairs / 1e6, protein_codinggenes, label = id)) +
  geom_point(colour = "steelblue", size = 2) + geom_smooth(method = "lm", se = FALSE, colour = "grey40") +
  geom_text(vjust = -0.7, size = 3) +
  labs(x = "Length (Mb)", y = "Protein-coding genes") + theme_minimal()
p_mi <- ggplot(chromosome, aes(basepairs / 1e6, mi_rna, label = id)) +
  geom_point(colour = "firebrick", size = 2) + geom_smooth(method = "lm", se = FALSE, colour = "grey40") +
  geom_text(vjust = -0.7, size = 3) +
  labs(x = "Length (Mb)", y = "miRNAs") + theme_minimal()
p_pc; p_mi
# Pearson correlations
cor(chromosome$basepairs, chromosome$protein_codinggenes)
cor(chromosome$basepairs, chromosome$mi_rna)
```

Both counts correlate positively with length (Pearson r ≈ `r round(cor(chromosome$basepairs, chromosome$protein_codinggenes), 2)` for protein-coding genes, ≈ `r round(cor(chromosome$basepairs, chromosome$mi_rna), 2)` for miRNAs), but with clear outliers: chr19 and chr17 are gene-dense, chr13, chr18 and chr X/Y are gene-poor for their size.

```{r}
# d. summary statistics for protein length and mass
proteins %>%
  summarise(across(c(length, mass), list(mean = mean, median = median, max = max), .names = "{.col}_{.fn}")) %>%
  pivot_longer(everything(), names_to = c("variable", "statistic"), names_pattern = "(.*)_(mean|median|max)") %>%
  pivot_wider(names_from = statistic, values_from = value)

# length vs mass on log-log axes (both right-skewed), with linear fit
ggplot(proteins, aes(length, mass / 1000)) +
  geom_point(alpha = 0.3, size = 0.8, colour = "#2c7fb8") +
  geom_smooth(method = "lm", colour = "firebrick", se = FALSE, linewidth = 0.6) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Protein length (amino acids, log10)", y = "Mass (kDa, log10)",
       title = "Protein mass scales linearly with length",
       subtitle = "20,430 human proteins; slope ≈ 0.11 kDa per residue") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
```

Mass is an almost deterministic linear function of length (~110 Da per residue); the log–log plot shows a straight line with few deviations (heavily modified or unusual-composition proteins). The distribution is right-skewed: median length 414 aa versus a maximum of 34,350 aa (titin).

# Session info

```{r}
sessionInfo()
```