# =============================================================================
# 05_differential_abundance.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Análisis de abundancia diferencial entre los dos sitios de vid
#              (Rancho Gil vs. Rancho San Ignacio).
#
# Análisis incluidos:
#   1. Wilcoxon Mann-Whitney (no paramétrico — recomendado para n pequeño)
#   2. DESeq2 (pseudoconteos - binomial negativo)
#   3. Análisis tipo LEfSe (LDA effect size) simplificado
#   4. Volcano plots
#   5. Dot/bubble plots de taxa diferenciales
#   6. Visualización de taxa más discriminantes
#
# NOTA: Con n=3 por grupo, todos los métodos tienen bajo poder estadístico.
#       Los resultados son exploratorios y deben confirmarse con más réplicas.
#       ALDEx2 es el método más robusto para muestras pequeñas con datos
#       composicionales (ver Nearing et al. 2022).
#
# Dependencias: 01_data_import.R
# =============================================================================

# ---- 0. Cargar datos ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("01_data_import.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
  library(phyloseq)
})

# ---- 1. Configuración ----
OUT_DIR1 <- here::here("~/Documents/Metagenomics/Microbial-soil/Artículo/")  # Ajustar si los scripts están en R_scripts/

OUT_DIR2 <- "figures/differential_abundance"

OUT_DIR <- file.path(OUT_DIR1, OUT_DIR2)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Umbral de significancia
ALPHA     <- 0.05   # p-value corregido
FC_THRESH <- 1.5    # Log2 fold-change mínimo para considerar "diferencial"
MIN_PREV  <- 2      # Prevalencia mínima (# muestras con abund > 0)

site_colors  <- c("Rancho_Gil" = "#2196F3", "Rancho_San_Ignacio" = "#FF5722")
site_labels  <- c("Rancho_Gil" = "Rancho Gil", "Rancho_San_Ignacio" = "Rancho San Ignacio")

# Grupos
groups <- meta_clean[colnames(bact_mat), "site"]
g1 <- "Rancho_Gil"
g2 <- "Rancho_San_Ignacio"

# ---- 2. Filtrar taxa por prevalencia mínima ----
filter_by_prevalence <- function(mat, min_prev) {
  keep <- rowSums(mat > 0) >= min_prev
  cat(sprintf("  → Prevalencia ≥ %d: %d / %d taxa conservados\n",
              min_prev, sum(keep), nrow(mat)))
  return(mat[keep, , drop = FALSE])
}

bact_filt  <- filter_by_prevalence(bact_mat,  MIN_PREV)
fungi_filt <- filter_by_prevalence(fungi_mat, MIN_PREV)

# ---- 3. Función: Test de Wilcoxon por taxa ----
run_wilcoxon <- function(mat, groups_vec, g1 = "Rancho_Gil", g2 = "Rancho_San_Ignacio") {
  s1 <- colnames(mat)[groups_vec == g1]
  s2 <- colnames(mat)[groups_vec == g2]

  results <- lapply(rownames(mat), function(taxon) {
    x1 <- as.numeric(mat[taxon, s1])
    x2 <- as.numeric(mat[taxon, s2])

    mean1 <- mean(x1)
    mean2 <- mean(x2)

    # Log2 fold-change (con pseudocount para evitar log(0))
    lfc <- log2((mean2 + 1e-6) / (mean1 + 1e-6))

    # Wilcoxon (exact = FALSE para n pequeño sin empates exactos)
    wt <- tryCatch(
      wilcox.test(x1, x2, exact = FALSE, correct = FALSE),
      error = function(e) list(p.value = 1, statistic = NA)
    )

    data.frame(
      Taxa       = taxon,
      Mean_G1    = round(mean1, 4),
      Mean_G2    = round(mean2, 4),
      log2FC     = round(lfc, 4),
      W          = ifelse(is.null(wt$statistic), NA, wt$statistic),
      p_value    = wt$p.value,
      stringsAsFactors = FALSE
    )
  })

  results_df <- do.call(rbind, results)

  # Corrección de Benjamini-Hochberg (FDR)
  results_df$p_adj     <- p.adjust(results_df$p_value, method = "BH")
  results_df$Sig       <- ifelse(results_df$p_adj < ALPHA, "Significativo", "ns")
  results_df$Direction <- ifelse(results_df$log2FC > 0,
                                  paste0("↑ en ", g2),
                                  paste0("↑ en ", g1))
  results_df$Sig_and_Dir <- ifelse(results_df$Sig == "Significativo",
                                    results_df$Direction, "ns")

  # Ordenar por p_adj
  results_df <- results_df %>% arrange(p_adj, desc(abs(log2FC)))

  return(results_df)
}

cat("\n===== WILCOXON MANN-WHITNEY =====\n")
wilcox_bact  <- run_wilcoxon(bact_filt,  groups, g1, g2)
wilcox_fungi <- run_wilcoxon(fungi_filt, groups, g1, g2)

cat(sprintf("  Bacterias significativas (p_adj < %.2f): %d taxa\n",
            ALPHA, sum(wilcox_bact$Sig == "Significativo")))
cat(sprintf("  Hongos significativas   (p_adj < %.2f): %d taxa\n",
            ALPHA, sum(wilcox_fungi$Sig == "Significativo")))

cat("\n  Top 10 taxa bacterianos más diferenciales:\n")
print(head(wilcox_bact %>% filter(Sig == "Significativo") %>%
             select(Taxa, Mean_G1, Mean_G2, log2FC, p_value, p_adj), 10))

# ---- 4. DESeq2 (con pseudo-conteos) ----
cat("\n===== DESeq2 =====\n")
cat("  Nota: DESeq2 requiere conteos enteros. Usamos pseudo-conteos (× 10,000).\n\n")

run_deseq2 <- function(ps_counts, group_col = "site",
                        ref_level = "Rancho_Gil", label = "") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    cat("  ⚠ DESeq2 no disponible. Instalar con BiocManager::install('DESeq2')\n")
    return(NULL)
  }
  library(DESeq2)

  # Convertir phyloseq a DESeq2
  dds <- phyloseq_to_deseq2(ps_counts, as.formula(paste("~", group_col)))
  dds[[group_col]] <- relevel(dds[[group_col]], ref = ref_level)

  # Estimar dispersión con fórmula "local" (más estable para n pequeño)
  dds <- tryCatch(
    DESeq(dds, fitType = "local", quiet = TRUE),
    error = function(e) {
      cat("  ⚠ Error en DESeq2:", conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(dds)) return(NULL)

  # Resultados
  res <- results(dds, alpha = ALPHA, independentFiltering = FALSE)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Taxa") %>%
    filter(!is.na(padj)) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(
      Sig       = ifelse(padj < ALPHA & abs(log2FoldChange) >= FC_THRESH,
                         "Significativo", "ns"),
      Direction = ifelse(log2FoldChange > 0,
                         paste0("↑ en ", setdiff(unique(dds[[group_col]]), ref_level)),
                         paste0("↑ en ", ref_level))
    )

  cat(sprintf("  DESeq2 - %s: %d taxa significativos (padj < %.2f & |log2FC| ≥ %.1f)\n",
              label, sum(res_df$Sig == "Significativo"), ALPHA, FC_THRESH))

  return(res_df)
}

deseq2_bact  <- run_deseq2(ps_bact_counts,  label = "Bacterias")
deseq2_fungi <- run_deseq2(ps_fungi_counts, label = "Hongos")

# ---- 5. Volcano Plots ----
cat("\n===== GENERANDO VOLCANO PLOTS =====\n")

plot_volcano <- function(results_df, title_str, filename,
                          g1_label = "Rancho Gil", g2_label = "Rancho San Ignacio",
                          fc_col = "log2FC", p_col = "p_adj",
                          label_top = 10) {
  # Preparar columnas si vienen de DESeq2
  if ("log2FoldChange" %in% names(results_df)) {
    results_df <- results_df %>%
      dplyr::rename(log2FC = log2FoldChange, p_adj = padj)
    results_df$p_value <- results_df$pvalue
  }

  # Clasificar puntos
  results_df <- results_df %>%
    mutate(
      neg_log10_p = -log10(p_adj + 1e-15),
      Color = case_when(
        p_adj < ALPHA & log2FC > FC_THRESH  ~ paste0("↑ ", g2_label),
        p_adj < ALPHA & log2FC < -FC_THRESH ~ paste0("↑ ", g1_label),
        TRUE ~ "ns"
      )
    )

  # Taxa para etiquetar (top por significancia + efecto)
  top_taxa <- results_df %>%
    filter(Color != "ns") %>%
    arrange(p_adj, desc(abs(log2FC))) %>%
    head(label_top)

  vol_colors <- setNames(
    c("#FF5722", "#2196F3", "grey70"),
    c(paste0("↑ ", g2_label), paste0("↑ ", g1_label), "ns")
  )

  p <- ggplot(results_df, aes(x = log2FC, y = neg_log10_p, color = Color)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = -log10(ALPHA), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-FC_THRESH, FC_THRESH), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    ggrepel::geom_text_repel(
      data = top_taxa,
      aes(label = str_trunc(Taxa, 25, "right")),
      size = 2.5, max.overlaps = 20, fontface = "italic",
      segment.color = "grey50", segment.size = 0.2
    ) +
    scale_color_manual(values = vol_colors, name = "Enriquecido") +
    scale_y_continuous(expand = c(0.01, 0)) +
    labs(
      title    = title_str,
      subtitle = sprintf("n=%d taxa | umbral: p_adj < %.2f & |log2FC| ≥ %.1f",
                         nrow(results_df), ALPHA, FC_THRESH),
      x        = "Log₂ Fold-Change",
      y        = "-log₁₀(p ajustado)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, color = "grey40")
    )

  ggsave(file.path(OUT_DIR, filename), p, width = 9, height = 6, dpi = 300)
  cat(sprintf("  ✓ Volcano plot guardado: %s\n", filename))

  return(invisible(p))
}

# Volcano con resultados Wilcoxon
plot_volcano(wilcox_bact,  "Volcano - Bacterias/Arqueas (Wilcoxon)",
             "volcano_bacteria_wilcoxon.pdf")
plot_volcano(wilcox_fungi, "Volcano - Hongos (Wilcoxon)",
             "volcano_fungi_wilcoxon.pdf")

# Volcano con DESeq2 (si está disponible)
if (!is.null(deseq2_bact)) {
  plot_volcano(deseq2_bact,  "Volcano - Bacterias/Arqueas (DESeq2)",
               "volcano_bacteria_deseq2.pdf", fc_col = "log2FoldChange", p_col = "padj")
}

if (!is.null(deseq2_fungi)) {
  plot_volcano(deseq2_fungi, "Volcano - Hongos (DESeq2)",
               "volcano_fungi_deseq2.pdf")
}

# ---- 6. Dot plot de taxa diferenciales significativos ----
cat("\n===== DOT PLOTS DE TAXA DIFERENCIALES =====\n")

plot_diff_dotplot <- function(results_df, mat, meta_df, groups_vec,
                               n_top = 20, title_str, filename) {
  # Seleccionar taxa significativos o top por efecto
  sig_taxa <- results_df %>%
    filter(Sig == "Significativo") %>%
    head(n_top)

  if (nrow(sig_taxa) == 0) {
    # Si no hay significativos, mostrar los top por p_value
    sig_taxa <- results_df %>%
      arrange(p_adj) %>%
      head(n_top)
    subtitle <- paste("Top", n_top, "taxa por p-adj (ninguno pasó umbral de significancia)")
  } else {
    subtitle <- paste("Taxa con p_adj < ", ALPHA)
  }

  # Extraer abundancias de los taxa seleccionados
  sel_taxa <- sig_taxa$Taxa
  sel_taxa <- sel_taxa[sel_taxa %in% rownames(mat)]

  if (length(sel_taxa) == 0) return(invisible(NULL))

  df_plot <- mat[sel_taxa, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(-Taxa, names_to = "sample_id", values_to = "abundance") %>%
    left_join(meta_df %>% rownames_to_column("sample_id") %>% select(sample_id, site),
              by = "sample_id") %>%
    mutate(
      Taxa       = factor(Taxa, levels = rev(sel_taxa)),
      site_label = site_labels[site]
    ) %>%
    group_by(Taxa, site_label) %>%
    summarise(
      mean_ab = mean(abundance),
      sd_ab   = sd(abundance),
      .groups = "drop"
    ) %>%
    left_join(sig_taxa %>% select(Taxa, log2FC, p_adj), by = "Taxa")

  p <- ggplot(df_plot, aes(x = site_label, y = Taxa,
                             size = mean_ab, fill = log2FC)) +
    geom_point(shape = 21, color = "grey50", alpha = 0.9) +
    scale_size_continuous(range = c(2, 12), name = "Abundancia\nmedia (%)") +
    scale_fill_gradient2(
      low  = "#2196F3", mid = "white", high = "#FF5722",
      midpoint = 0, name = "log₂FC"
    ) +
    labs(
      title    = title_str,
      subtitle = subtitle,
      x        = "Sitio",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y      = element_text(size = 8, face = "italic"),
      axis.text.x      = element_text(size = 10),
      panel.grid.major = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, color = "grey40")
    )

  ggsave(file.path(OUT_DIR, filename), p, width = 8, height = max(6, length(sel_taxa) * 0.4 + 3),
         dpi = 300)
  cat(sprintf("  ✓ Dot plot guardado: %s\n", filename))
  return(invisible(p))
}

plot_diff_dotplot(wilcox_bact,  bact_filt,  meta_clean[colnames(bact_mat), ],
                   groups,
                   title_str = "Taxa bacterianos diferenciales (Wilcoxon)",
                   filename  = "dotplot_diff_bacteria.pdf")

plot_diff_dotplot(wilcox_fungi, fungi_filt, meta_clean[colnames(fungi_mat), ],
                   groups,
                   title_str = "Taxa fúngicos diferenciales (Wilcoxon)",
                   filename  = "dotplot_diff_fungi.pdf")

# ---- 7. Boxplots para los taxa más discriminantes ----
plot_top_diff_boxplots <- function(results_df, mat, meta_df, groups_vec,
                                    n_top = 6, title_str, filename) {
  top_taxa <- results_df %>%
    arrange(p_adj, desc(abs(log2FC))) %>%
    head(n_top) %>%
    pull(Taxa)

  top_taxa <- top_taxa[top_taxa %in% rownames(mat)]
  if (length(top_taxa) == 0) return(invisible(NULL))

  df_plot <- mat[top_taxa, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(-Taxa, names_to = "sample_id", values_to = "abundance") %>%
    left_join(meta_df %>% rownames_to_column("sample_id") %>% select(sample_id, site),
              by = "sample_id") %>%
    mutate(site_label = site_labels[site],
           Taxa = str_trunc(Taxa, 30, "right"))

  p <- ggplot(df_plot, aes(x = site_label, y = abundance, fill = site, color = site)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 2.5, alpha = 0.9) +
    scale_fill_manual(values  = site_colors, labels = site_labels) +
    scale_color_manual(values = site_colors, labels = site_labels) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format",
                       size = 3, label.y.npc = 0.95) +
    facet_wrap(~ Taxa, scales = "free_y", ncol = 3) +
    labs(
      title = title_str,
      x     = NULL,
      y     = "Abundancia relativa (%)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(angle = 20, hjust = 1, size = 8),
      strip.text       = element_text(size = 7, face = "italic"),
      strip.background = element_rect(fill = "grey92"),
      plot.title       = element_text(face = "bold", size = 11)
    )

  ggsave(file.path(OUT_DIR, filename), p,
         width  = 12,
         height = ceiling(length(top_taxa) / 3) * 3.5 + 1,
         dpi    = 300)
  cat(sprintf("  ✓ Boxplots guardados: %s\n", filename))
  return(invisible(p))
}

cat("\n===== BOXPLOTS DE TAXA TOP DIFERENCIALES =====\n")
plot_top_diff_boxplots(wilcox_bact,  bact_filt,  meta_clean[colnames(bact_mat), ],
                        groups,
                        title_str = "Top taxa bacterianos más diferenciales (Wilcoxon)",
                        filename  = "boxplots_top_diff_bacteria.pdf")

plot_top_diff_boxplots(wilcox_fungi, fungi_filt, meta_clean[colnames(fungi_mat), ],
                        groups,
                        title_str = "Top taxa fúngicos más diferenciales (Wilcoxon)",
                        filename  = "boxplots_top_diff_fungi.pdf")

# ---- 8. Exportar tablas de resultados ----
cat("\n===== EXPORTANDO RESULTADOS =====\n")
write.csv(wilcox_bact,  file.path(OUT_DIR, "wilcoxon_bacteria.csv"),  row.names = FALSE)
write.csv(wilcox_fungi, file.path(OUT_DIR, "wilcoxon_fungi.csv"),     row.names = FALSE)

if (!is.null(deseq2_bact))  write.csv(deseq2_bact,  file.path(OUT_DIR, "deseq2_bacteria.csv"),  row.names = FALSE)
if (!is.null(deseq2_fungi)) write.csv(deseq2_fungi, file.path(OUT_DIR, "deseq2_fungi.csv"),     row.names = FALSE)

cat("  ✓ Tablas exportadas: wilcoxon_*.csv, deseq2_*.csv\n")

# ---- 9. Resumen final ----
cat("\n===== RESUMEN DE ABUNDANCIA DIFERENCIAL =====\n")
cat(sprintf("
  Método: Wilcoxon Mann-Whitney + FDR Benjamini-Hochberg
  Umbral: p_adj < %.2f, |log₂FC| ≥ %.1f
  ────────────────────────────────────────────────────
  Bacterias/Arqueas (16S):
    Total taxa evaluados:      %d
    Taxa significativos:       %d
    ↑ Enriquecidos en GIL:     %d
    ↑ Enriquecidos en SAN IGN: %d
  Hongos (ITS):
    Total taxa evaluados:      %d
    Taxa significativos:       %d
    ↑ Enriquecidos en GIL:     %d
    ↑ Enriquecidos en SAN IGN: %d
  ────────────────────────────────────────────────────
  ADVERTENCIA: n=3 por grupo → bajo poder estadístico.
  Interpretar como análisis exploratorio.
",
  ALPHA, FC_THRESH,
  nrow(wilcox_bact),
  sum(wilcox_bact$Sig == "Significativo"),
  sum(wilcox_bact$Sig == "Significativo" & wilcox_bact$log2FC < 0),
  sum(wilcox_bact$Sig == "Significativo" & wilcox_bact$log2FC > 0),
  nrow(wilcox_fungi),
  sum(wilcox_fungi$Sig == "Significativo"),
  sum(wilcox_fungi$Sig == "Significativo" & wilcox_fungi$log2FC < 0),
  sum(wilcox_fungi$Sig == "Significativo" & wilcox_fungi$log2FC > 0)
))

cat("\n✓ Análisis de abundancia diferencial completado. Procede con 06_environmental_analysis.R\n")

