# Background

In this script we are going to examine the richness of biosynthetic gene
clusters (BGCs) within some cyanobacterial genomes housed within the
BLCC. Genomes have already been assembled and annotated. Here we are
going to run them through
[antiSMASH](https://antismash.secondarymetabolites.org/#!/about), using
[multiSMASH](https://github.com/zreitz/multismash), and explore those
data. Additionally, we are going to build a quick phylogenetic tree with
[IQ-TREE](http://www.iqtree.org/) using genes identified with
[GToTree](https://github.com/AstrobioMike/GToTree).

# Run multiSMASH

This isn’t run in R, I ran this on a cluster.

``` r
conda activate antismash7
N=MULTISMASH
CMD="multismash config.yaml"
sbatch -A NAME -J ${N} -c 64 --mem=500G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=ALL --mail-user=ME -t 48:00:00 --wrap="${CMD}"
```

# Make phylogenetic tree

Again, this isn’t run in R, I ran this on a cluster.

## Identify genes

``` r
conda activate gtotree
N=GTOTREE
CMD="GToTree -f genomelist.txt -H Cyanobacteria -B -D -n 2 -j 16 -M 2 -N"
sbatch -A NAME -J ${N} -c 32 --mem=200G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=ME -t 48:00:00 --wrap="${CMD}"
```

# Construct the phylogeny

``` r
N=GTOTREE
CMD="iqtree -s GToTree_output/Aligned_SCGs.faa -T 64 -m LG+I+G4 -B 1000 --seed 42069 -pre iqtree_out2"
sbatch -A NAME -J ${N} -c 64 --mem=200G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=ME -t 48:00:00 --wrap="${CMD}"
```

# Now were in R!

# Load Packages

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggalign)
```

    ## Warning: package 'ggalign' was built under R version 4.3.3

``` r
library(RColorBrewer)
library(ape)
```

    ## 
    ## Attaching package: 'ape'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

``` r
library(ggtree)
```

    ## ggtree v3.10.1 For help: https://yulab-smu.top/treedata-book/
    ## 
    ## If you use the ggtree package suite in published research, please cite
    ## the appropriate paper(s):
    ## 
    ## Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
    ## ggtree: an R package for visualization and annotation of phylogenetic
    ## trees with their covariates and other associated data. Methods in
    ## Ecology and Evolution. 2017, 8(1):28-36. doi:10.1111/2041-210X.12628
    ## 
    ## Guangchuang Yu. Using ggtree to visualize data on tree-like structures.
    ## Current Protocols in Bioinformatics. 2020, 69:e96. doi:10.1002/cpbi.96
    ## 
    ## Shuangbin Xu, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan
    ## Dai, Tommy T. Lam, Yi Guan, Guangchuang Yu. Ggtree: A serialized data
    ## object for visualization of a phylogenetic tree and annotation data.
    ## iMeta 2022, 1(4):e56. doi:10.1002/imt2.56
    ## 
    ## Attaching package: 'ggtree'
    ## 
    ## The following object is masked from 'package:ape':
    ## 
    ##     rotate
    ## 
    ## The following object is masked from 'package:ggalign':
    ## 
    ##     inset
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

## Make a quick function

Function to get the first word before a space, dash, or slash

``` r
get_first_part = function(x) {
  word = str_extract(x, "^[^\\s/-]+")
  return(ifelse(is.na(word), x, word))
}
```

# Read in the data

# Plot the tree

## Leaf order from tree

``` r
leaf_order = tree$tip.label
```

### Another leave order

``` r
# Create a named vector for easy lookup
name_lookup <- setNames(names$new_name, names$old_name)

# Create leaf_order2 by replacing values
leaf_order2 <- sapply(leaf_order, function(x) {
  if (x %in% names(name_lookup)) {
    return(name_lookup[x])
  } else {
    return(x)
  }
})

# If you want leaf_order2 to be a simple vector instead of a named vector
leaf_order2 <- unname(leaf_order2)
```

## Plot the thing

``` r
tree_plot <- ggtree(tree) + 
  geom_nodelab(aes(label=label, x=branch), size=3, vjust=-0.5) +
  theme_tree2() +
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
```

# Massage the data

Some of the names are way to long, that really annoying for plotting so
were going to trim them down

``` r
# Extract the unique file names from the BGCs data
bgcs_files <- unique(df$file)

# Find the common entries
common_entries <- intersect(leaf_order, bgcs_files)

# Filter the BGCs data to include only the common entries
filtered_bgcs_data <- df %>%
  filter(file %in% common_entries)

transformed_df = filtered_bgcs_data %>%
  #filter(!is.na(KCB_hit) & KCB_hit != "") %>%
  mutate(
    KCB_hit = str_extract(KCB_hit, "^[^/]+"),  # Extract everything before the first "/"
    value = 1
  ) %>%
  mutate(KCB_hit = str_replace(KCB_hit, "anabaenopeptin (788|908|NZ857).*", "anabaenopeptin"))%>%
  mutate(KCB_hit = str_replace(KCB_hit, "anacyclamide (A10|D8P).*", "anacyclamide"))%>%
  mutate(KCB_hit = str_replace(KCB_hit, "malyngamide (C acetate|I).*", "malyngamide"))%>%
  mutate(KCB_hit = str_replace(KCB_hit, "microcystin LR", "microcystin"))%>%
  mutate(KCB_hit = str_replace(KCB_hit, "microcystin-LR", "microcystin"))%>%
  mutate(KCB_hit = str_replace(KCB_hit, "microviridin (1688|1777|B|K|N9).*", "microviridin"))%>%
  mutate(KCB_hit = str_replace(KCB_hit, "nostopeptolide (A1|A2).*", "nostopeptolide"))%>%
  mutate(KCB_hit = str_replace(KCB_hit, "puwainaphycin (A|F).*", "puwainaphycin"))
```

## BGCs of interest

These are the ones of interest. Why? because I said so. Some are
cyanotoxins, some have cool names. I prob missed some cool ones, another
problem for another day.

``` r
kcb_hit_list = c(
  "aeruginoside 126B", "aeruginosin 98-A", "amycomicin", "anabaenopeptin", "anatoxin-a",
  "apratoxin A", "anachelin", "cyanochelin A", "cyanopeptolin", "cylindrocyclophane",
  "cylindrocyclophane D", "cylindrospermopsin", "cryptophycin-327", "cyphomycin",
  "curacin A", "microcystin", "geosmin", "jamaicamide A", "lankacidin C",
  "malyngamide", "microginin", "micropeptin K139", "microviridin", "nocuolin A",
  "nostoclide N1", "nostocyclopeptide A2", "nostolysamide A",
  "nostopeptolide", "nostophycin", "oxalomycin B", "pseudospumigin A",
  "puwainaphycin", "saxitoxin")
```

# Heat map

## Prepare the data for heatmap

``` r
heatmap_data = transformed_df %>%
  mutate(KCB_hit = str_extract(KCB_hit, "[^/]+")) %>%  # Extract everything before the first "/"
  mutate(KCB_hit = ifelse(KCB_hit %in% kcb_hit_list, KCB_hit, "Other")) %>%
  distinct(file, KCB_hit) %>%
  mutate(value = 1) %>%
  pivot_wider(
    id_cols = file,
    names_from = KCB_hit,
    values_from = value,
    values_fill = list(value = 0)
  ) %>%
  pivot_longer(cols = -file, names_to = "KCB_hit", values_to = "value") %>%
  mutate(file_type = case_when(
    str_starts(file, "F") ~ "Freshwater",
    str_starts(file, "M") ~ "Marine",
    str_starts(file, "T") ~ "Terrestrial",
    TRUE ~ "Other"
  )) 
```

## Replace only the matching values in the ‘file’ column

Lets get some better names, shall we?

``` r
heatmap_data2 = heatmap_data %>%
  mutate(file = factor(file, levels = leaf_order)) %>%
  arrange(file) %>% #order by tree order
  left_join(names, by = c("file" = "old_name")) %>%
  mutate(file = coalesce(new_name, file)) %>%
  select(-new_name)  # Remove the new_name column after using it
```

## Create the heatmap

``` r
heatmap_plot = ggplot(heatmap_data2, aes(x = KCB_hit, y = file)) +
  geom_tile(aes(fill = factor(value)), color = "black", width = 1) +
  geom_tile(aes(x = -1, fill = file_type), width = 2) +  # New column for file type
  scale_fill_manual(values = c("0" = "white", "1" = "black",
                               "Freshwater" = "green", "Marine" = "blue", 
                               "Terrestrial" = "brown", "Other" = "grey"),
                    name = "Legend",
                    labels = c("0" = "Absent", "1" = "Present",
                      "Freshwater" = "Freshwater", 
                      "Marine" = "Marine", "Terrestrial" = "Terrestrial", 
                      "Other" = "Other files")
                    ) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +  # Increase legend key size
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8)
  ) +
  scale_x_discrete(expand = c(0, 0)) +  # Remove padding on x-axis
  scale_y_discrete(limits = leaf_order2, expand = c(0, 0))    # Remove padding on y-axis
```

# Bar plot

Now lets make a stacked bar plot of the number of BGCs for each genome

## Prepare data for stacked bar plot with collapsed product categories

There are quite a few differnt types of BGCs, were going to pick the 15
most common ones and look at them

``` r
bar_data = df %>%
  filter(!is.na(KCB_hit) & KCB_hit != "") %>%
  #mutate(product = get_first_part(product)) %>%
  group_by(file, product) %>%
  summarize(count = n(), .groups = "drop") %>%
  mutate(product = fct_lump(product, n = 15, w = count))  # Group less common products as "Other"
```

## Replace the names, again

``` r
bar_data2 = bar_data %>%
  left_join(names, by = c("file" = "old_name")) %>%
  mutate(file = coalesce(new_name, file)) %>%
  select(-new_name)  # Remove the new_name column after using it
```

## Define a color palette with 16 distinct colors

``` r
color_palette = c(
  brewer.pal(7, "Set2"),
  brewer.pal(9, "Set1")
)
```

## Create the stacked bar plot with legend on the right and new color palette

``` r
bar_plot = ggplot(bar_data2, aes(x = count, y = file, fill = product)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() + 
  scale_y_discrete(limits = leaf_order2, expand = c(0, 0), position = "right") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank()
        ) +
  guides(fill = guide_legend(nrow = 2))
```

# How many unclassified BGCs?

Theres still some novelty out there. These BGCs were not classified as
anything within [MiBiG](https://mibig.secondarymetabolites.org/)

## Prepare data for NA count

``` r
na_count_data = df %>%
  group_by(file) %>%
  summarize(na_count = sum(is.na(KCB_hit) | KCB_hit == ""))
```

## Replace only the matching values in the ‘file’ column

``` r
na_count_data2 = na_count_data %>%
  left_join(names, by = c("file" = "old_name")) %>%
  mutate(file = coalesce(new_name, file)) %>%
  select(-new_name)  # Remove the new_name column after using it
```

## Create the NA count plot as a single column

``` r
na_count_plot = ggplot(na_count_data2, aes(y = file)) +
  geom_text(aes(label = na_count, x = 0), size = 3) + # Adjust size here
  theme_void() +
  theme(plot.margin = margin(5.5, 0, 5.5, 0))+ 
  scale_y_discrete(limits = leaf_order2, expand = c(0, 0))
```

# Lets tie it all together now

Align the plots using ggalign

``` r
aligned_plot <- align_plots(
  tree_plot,
  heatmap_plot,
  na_count_plot,
  bar_plot,
  widths = c(0.8,2, 0.1, 1)
)
```

# Session info

``` r
devtools::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.3.2 (2023-10-31)
    ##  os       macOS 15.0.1
    ##  system   aarch64, darwin20
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       America/New_York
    ##  date     2024-10-25
    ##  pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package      * version date (UTC) lib source
    ##  ape          * 5.8     2024-04-11 [1] CRAN (R 4.3.1)
    ##  aplot          0.2.3   2024-06-17 [1] CRAN (R 4.3.3)
    ##  bit            4.5.0   2024-09-20 [1] CRAN (R 4.3.3)
    ##  bit64          4.5.2   2024-09-22 [1] CRAN (R 4.3.3)
    ##  cachem         1.1.0   2024-05-16 [1] CRAN (R 4.3.3)
    ##  cli            3.6.3   2024-06-21 [1] CRAN (R 4.3.3)
    ##  colorspace     2.1-1   2024-07-26 [1] CRAN (R 4.3.3)
    ##  crayon         1.5.3   2024-06-20 [1] CRAN (R 4.3.3)
    ##  data.table     1.16.2  2024-10-10 [1] CRAN (R 4.3.3)
    ##  devtools       2.4.5   2022-10-11 [1] CRAN (R 4.3.0)
    ##  digest         0.6.37  2024-08-19 [1] CRAN (R 4.3.3)
    ##  dplyr        * 1.1.4   2023-11-17 [1] CRAN (R 4.3.1)
    ##  ellipsis       0.3.2   2021-04-29 [1] CRAN (R 4.3.0)
    ##  evaluate       1.0.1   2024-10-10 [1] CRAN (R 4.3.3)
    ##  fansi          1.0.6   2023-12-08 [1] CRAN (R 4.3.1)
    ##  farver         2.1.2   2024-05-13 [1] CRAN (R 4.3.3)
    ##  fastmap        1.2.0   2024-05-15 [1] CRAN (R 4.3.3)
    ##  forcats      * 1.0.0   2023-01-29 [1] CRAN (R 4.3.0)
    ##  fs             1.6.4   2024-04-25 [1] CRAN (R 4.3.1)
    ##  generics       0.1.3   2022-07-05 [1] CRAN (R 4.3.0)
    ##  ggalign      * 0.0.4   2024-10-12 [1] CRAN (R 4.3.3)
    ##  ggfun          0.1.6   2024-08-28 [1] CRAN (R 4.3.3)
    ##  ggplot2      * 3.5.1   2024-04-23 [1] CRAN (R 4.3.1)
    ##  ggplotify      0.1.2   2023-08-09 [1] CRAN (R 4.3.0)
    ##  ggtree       * 3.10.1  2024-02-27 [1] Bioconductor 3.18 (R 4.3.2)
    ##  glue           1.8.0   2024-09-30 [1] CRAN (R 4.3.3)
    ##  gridGraphics   0.5-1   2020-12-13 [1] CRAN (R 4.3.0)
    ##  gtable         0.3.5   2024-04-22 [1] CRAN (R 4.3.1)
    ##  hms            1.1.3   2023-03-21 [1] CRAN (R 4.3.0)
    ##  htmltools      0.5.8.1 2024-04-04 [1] CRAN (R 4.3.1)
    ##  htmlwidgets    1.6.4   2023-12-06 [1] CRAN (R 4.3.1)
    ##  httpuv         1.6.15  2024-03-26 [1] CRAN (R 4.3.1)
    ##  jsonlite       1.8.9   2024-09-20 [1] CRAN (R 4.3.3)
    ##  knitr          1.48    2024-07-07 [1] CRAN (R 4.3.3)
    ##  later          1.3.2   2023-12-06 [1] CRAN (R 4.3.1)
    ##  lattice        0.22-6  2024-03-20 [1] CRAN (R 4.3.1)
    ##  lazyeval       0.2.2   2019-03-15 [1] CRAN (R 4.3.0)
    ##  lifecycle      1.0.4   2023-11-07 [1] CRAN (R 4.3.1)
    ##  lubridate    * 1.9.3   2023-09-27 [1] CRAN (R 4.3.1)
    ##  magrittr       2.0.3   2022-03-30 [1] CRAN (R 4.3.0)
    ##  memoise        2.0.1   2021-11-26 [1] CRAN (R 4.3.0)
    ##  mime           0.12    2021-09-28 [1] CRAN (R 4.3.0)
    ##  miniUI         0.1.1.1 2018-05-18 [1] CRAN (R 4.3.0)
    ##  munsell        0.5.1   2024-04-01 [1] CRAN (R 4.3.1)
    ##  nlme           3.1-166 2024-08-14 [1] CRAN (R 4.3.3)
    ##  patchwork      1.3.0   2024-09-16 [1] CRAN (R 4.3.3)
    ##  pillar         1.9.0   2023-03-22 [1] CRAN (R 4.3.0)
    ##  pkgbuild       1.4.4   2024-03-17 [1] CRAN (R 4.3.1)
    ##  pkgconfig      2.0.3   2019-09-22 [1] CRAN (R 4.3.0)
    ##  pkgload        1.4.0   2024-06-28 [1] CRAN (R 4.3.3)
    ##  profvis        0.4.0   2024-09-20 [1] CRAN (R 4.3.3)
    ##  promises       1.3.0   2024-04-05 [1] CRAN (R 4.3.1)
    ##  purrr        * 1.0.2   2023-08-10 [1] CRAN (R 4.3.0)
    ##  R6             2.5.1   2021-08-19 [1] CRAN (R 4.3.0)
    ##  RColorBrewer * 1.1-3   2022-04-03 [1] CRAN (R 4.3.0)
    ##  Rcpp           1.0.13  2024-07-17 [1] CRAN (R 4.3.3)
    ##  readr        * 2.1.5   2024-01-10 [1] CRAN (R 4.3.1)
    ##  remotes        2.5.0   2024-03-17 [1] CRAN (R 4.3.1)
    ##  rlang          1.1.4   2024-06-04 [1] CRAN (R 4.3.3)
    ##  rmarkdown      2.28    2024-08-17 [1] CRAN (R 4.3.3)
    ##  rstudioapi     0.16.0  2024-03-24 [1] CRAN (R 4.3.1)
    ##  scales         1.3.0   2023-11-28 [1] CRAN (R 4.3.1)
    ##  sessioninfo    1.2.2   2021-12-06 [1] CRAN (R 4.3.0)
    ##  shiny          1.9.1   2024-08-01 [1] CRAN (R 4.3.3)
    ##  stringi        1.8.4   2024-05-06 [1] CRAN (R 4.3.1)
    ##  stringr      * 1.5.1   2023-11-14 [1] CRAN (R 4.3.1)
    ##  tibble       * 3.2.1   2023-03-20 [1] CRAN (R 4.3.0)
    ##  tidyr        * 1.3.1   2024-01-24 [1] CRAN (R 4.3.1)
    ##  tidyselect     1.2.1   2024-03-11 [1] CRAN (R 4.3.1)
    ##  tidytree       0.4.6   2023-12-12 [1] CRAN (R 4.3.1)
    ##  tidyverse    * 2.0.0   2023-02-22 [1] CRAN (R 4.3.0)
    ##  timechange     0.3.0   2024-01-18 [1] CRAN (R 4.3.1)
    ##  treeio         1.26.0  2023-11-06 [1] Bioconductor
    ##  tzdb           0.4.0   2023-05-12 [1] CRAN (R 4.3.0)
    ##  urlchecker     1.0.1   2021-11-30 [1] CRAN (R 4.3.0)
    ##  usethis        3.0.0   2024-07-29 [1] CRAN (R 4.3.3)
    ##  utf8           1.2.4   2023-10-22 [1] CRAN (R 4.3.1)
    ##  vctrs          0.6.5   2023-12-01 [1] CRAN (R 4.3.1)
    ##  vroom          1.6.5   2023-12-05 [1] CRAN (R 4.3.1)
    ##  withr          3.0.1   2024-07-31 [1] CRAN (R 4.3.3)
    ##  xfun           0.48    2024-10-03 [1] CRAN (R 4.3.3)
    ##  xtable         1.8-4   2019-04-21 [1] CRAN (R 4.3.0)
    ##  yaml           2.3.10  2024-07-26 [1] CRAN (R 4.3.3)
    ##  yulab.utils    0.1.7   2024-08-26 [1] CRAN (R 4.3.3)
    ## 
    ##  [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
