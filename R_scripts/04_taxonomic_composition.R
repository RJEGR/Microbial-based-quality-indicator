# =============================================================================
# 04_taxonomic_composition.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Análisis y visualización de la composición taxonómica de las
#              comunidades microbianas del suelo de vid.
#
# Análisis incluidos:
#   1. Barplots de composición relativa (Top 20 taxa / géneros)
#   2. Heatmaps de abundancia (pheatmap)
#   3. Diagramas de Venn / UpSet (taxa compartidos entre sitios)
#   4. Composición por género (agrupación del primer nivel taxonómico)
#   5. Burbuja (bubble chart) de taxa más abundantes
#
# Dependencias: 01_data_import.R
# =============================================================================

# ---- 0. Cargar datos ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("01_data_import.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(forcats)
  library(reshape2)
  library(pheatmap)
  library(RColorBrewer)
  library(patchwork)
  library(stringr)
  library(ggVennDiagram)
  library(scales)
})

# ---- 1. Configuración ----
OUT_DIR1 <- here::here("~/Documents/Metagenomics/Microbial-soil/Artículo/")  # Ajustar si los scripts están en R_scripts/

OUT_DIR2 <- "figures/taxonomic_composition"

OUT_DIR <- file.path(OUT_DIR1, OUT_DIR2)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

site_colors  <- c("Rancho_Gil" = "#2196F3", "Rancho_San_Ignacio" = "#FF5722")
site_labels  <- c("Rancho_Gil" = "Rancho Gil", "Rancho_San_Ignacio" = "Rancho San Ignacio")
TOP_N_TAXA   <- 20   # Número de taxa para barplots y heatmaps
TOP_N_GENERA <- 15   # Número de géneros

# ---- 2. Función: extraer Top-N taxa y colapsar el resto como "Otros" ----
get_top_taxa <- function(mat, n = 20, group_var = NULL) {
  # mat: abundancias relativas (taxa en filas, muestras en columnas)

  # Calcular abundancia media por taxa
  mean_abund <- rowMeans(mat)

  # Seleccionar Top-N
  top_taxa <- names(sort(mean_abund, decreasing = TRUE))[1:min(n, nrow(mat))]

  # Separar top y "Otros"
  mat_top   <- mat[top_taxa, , drop = FALSE]
  mat_other <- mat[!rownames(mat) %in% top_taxa, , drop = FALSE]
  other_row  <- colSums(mat_other)

  mat_final <- rbind(mat_top, "Otros" = other_row)

  return(mat_final)
}

# ---- 3. Función: agrupar por género ----
# Extrae el primer componente del nombre de especie como género
aggregate_by_genus <- function(mat) {
  taxa_names <- rownames(mat)
  genera      <- str_extract(taxa_names, "^[A-Za-z]+")

  # Agrupar taxa del mismo género sumando abundancias
  mat_df <- as.data.frame(mat) %>%
    mutate(Genus = genera) %>%
    group_by(Genus) %>%
    summarise(across(everything(), sum), .groups = "drop") %>%
    column_to_rownames("Genus") %>%
    as.matrix()

  return(mat_df)
}

# ---- 4. Preparar datos de composición ----
# Bacterias - a nivel de especie
bact_top  <- get_top_taxa(bact_mat,  n = TOP_N_TAXA)
fungi_top <- get_top_taxa(fungi_mat, n = TOP_N_TAXA)

# Bacterias - agrupadas por género
bact_genus  <- aggregate_by_genus(bact_mat)
fungi_genus <- aggregate_by_genus(fungi_mat)
bact_genus_top  <- get_top_taxa(bact_genus,  n = TOP_N_GENERA)
fungi_genus_top <- get_top_taxa(fungi_genus, n = TOP_N_GENERA)

cat(sprintf("  ✓ Top %d taxa bacterias:  %d (+ Otros)\n", TOP_N_TAXA, TOP_N_TAXA))
cat(sprintf("  ✓ Top %d taxa hongos:     %d (+ Otros)\n", TOP_N_TAXA, TOP_N_TAXA))
cat(sprintf("  ✓ Géneros únicos bacterias: %d\n", nrow(bact_genus)))
cat(sprintf("  ✓ Géneros únicos hongos:    %d\n", nrow(fungi_genus)))

# ---- 5. Función: barplot de composición ----
plot_composition_barplot <- function(mat, meta_df, title_str,
                                     facet_by = "site", nrow_legend = 5) {
  # Convertir a formato largo
  df <- as.data.frame(mat) %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(-Taxa, names_to = "sample_id", values_to = "abundance") %>%
    left_join(meta_df %>% rownames_to_column("sample_id") %>%
                select(sample_id, site, site_name, replicate),
              by = "sample_id") %>%
    mutate(
      Taxa = factor(Taxa, levels = c(rev(rownames(mat)[-nrow(mat)]), "Otros")),
      site_label = site_labels[site]
    )

  # Paleta de colores para taxa
  n_taxa  <- nrow(mat)
  taxa_colors <- c(
    colorRampPalette(brewer.pal(12, "Set3"))(n_taxa - 1),
    "grey70"  # color para "Otros"
  )
  names(taxa_colors) <- c(rownames(mat)[-nrow(mat)], "Otros")

  p <- ggplot(df, aes(x = sample_id, y = abundance, fill = Taxa)) +
    geom_bar(stat = "identity", width = 0.85, color = "white", linewidth = 0.1) +
    scale_fill_manual(values = taxa_colors) +
    scale_y_continuous(expand = c(0, 0), labels = label_number(suffix = "%")) +
    facet_wrap(~ site_label, scales = "free_x", nrow = 1) +
    labs(
      title = title_str,
      x     = "Muestra",
      y     = "Abundancia relativa (%)",
      fill  = "Taxa"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x       = element_text(angle = 35, hjust = 1, size = 9),
      legend.position   = "right",
      legend.text       = element_text(size = 7, face = "italic"),
      legend.title      = element_text(size = 9, face = "bold"),
      legend.key.size   = unit(0.4, "cm"),
      panel.grid        = element_blank(),
      strip.background  = element_rect(fill = "grey92"),
      strip.text        = element_text(face = "bold", size = 10),
      plot.title        = element_text(face = "bold", size = 12)
    ) +
    guides(fill = guide_legend(ncol = 1, reverse = TRUE))

  return(p)
}

# ---- 6. Generar barplots ----
cat("\n===== GENERANDO BARPLOTS DE COMPOSICIÓN =====\n")

# Nivel de especie - Top 20
p_bact_sp <- plot_composition_barplot(
  mat       = bact_top,
  meta_df   = meta_clean,
  title_str = paste0("Composición bacteriana — Top ", TOP_N_TAXA, " taxa (nivel especie)")
)

p_fungi_sp <- plot_composition_barplot(
  mat       = fungi_top,
  meta_df   = meta_clean,
  title_str = paste0("Composición fúngica — Top ", TOP_N_TAXA, " taxa (nivel especie)")
)

# Nivel de género - Top 15
p_bact_gen <- plot_composition_barplot(
  mat       = bact_genus_top,
  meta_df   = meta_clean,
  title_str = paste0("Composición bacteriana — Top ", TOP_N_GENERA, " géneros")
)

p_fungi_gen <- plot_composition_barplot(
  mat       = fungi_genus_top,
  meta_df   = meta_clean,
  title_str = paste0("Composición fúngica — Top ", TOP_N_GENERA, " géneros")
)

# Guardar
for (fig_name in c("bact_sp", "bact_gen", "fungi_sp", "fungi_gen")) {
  fig <- get(paste0("p_", fig_name))
  ggsave(file.path(OUT_DIR, paste0("barplot_", fig_name, ".pdf")),
         fig, width = 12, height = 7, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0("barplot_", fig_name, ".png")),
         fig, width = 12, height = 7, dpi = 300)
}
cat("  ✓ Barplots guardados (especie y género, bacterias y hongos)\n")

# ---- 7. Heatmaps de abundancia (pheatmap) ----
cat("\n===== GENERANDO HEATMAPS DE ABUNDANCIA =====\n")

plot_abund_heatmap <- function(mat, meta_df, title_str, n_top = 30,
                                filename = "heatmap.pdf") {
  # Seleccionar top taxa por abundancia media
  top_taxa <- names(sort(rowMeans(mat), decreasing = TRUE))[1:min(n_top, nrow(mat))]
  mat_top  <- mat[top_taxa, , drop = FALSE]

  # Escalar por fila (z-score) para comparar perfiles entre muestras
  mat_scaled <- t(scale(t(mat_top)))
  mat_scaled[is.nan(mat_scaled)] <- 0

  # Anotación de columnas (muestras)
  annot_col <- meta_df %>%
    select(site) %>%
    rename(Sitio = site)
  annot_col$Sitio <- as.character(annot_col$Sitio)

  ann_colors <- list(
    Sitio = c("Rancho_Gil" = "#2196F3", "Rancho_San_Ignacio" = "#FF5722")
  )

  # Paleta de color para heatmap (divergente: azul → blanco → rojo)
  heatmap_colors <- colorRampPalette(
    rev(brewer.pal(11, "RdBu"))
  )(100)

  pdf(file.path(OUT_DIR, filename), width = 10, height = 12)
  pheatmap(
    mat_scaled,
    color             = heatmap_colors,
    annotation_col    = annot_col,
    annotation_colors = ann_colors,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    clustering_method = "ward.D2",
    scale             = "none",       # ya escalamos manualmente
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize_row      = 7,
    fontsize_col      = 9,
    fontface_row      = "italic",
    border_color      = NA,
    main              = paste(title_str, "\n(z-score por taxa)"),
    angle_col         = 45,
    treeheight_col    = 30,
    treeheight_row    = 30
  )
  dev.off()
  cat(sprintf("  ✓ Heatmap guardado: %s\n", filename))
}

plot_abund_heatmap(bact_mat,  meta_clean[colnames(bact_mat), ],
                   "Top 30 Bacterias/Arqueas (16S)",  n_top = 30,
                   filename = "heatmap_bacteria_top30.pdf")

plot_abund_heatmap(fungi_mat, meta_clean[colnames(fungi_mat), ],
                   "Top 30 Hongos (ITS)",             n_top = 30,
                   filename = "heatmap_fungi_top30.pdf")

plot_abund_heatmap(bact_genus,  meta_clean[colnames(bact_mat), ],
                   "Top 30 Géneros Bacterianos (16S)", n_top = 30,
                   filename = "heatmap_bacteria_genera.pdf")

plot_abund_heatmap(fungi_genus, meta_clean[colnames(fungi_mat), ],
                   "Top 30 Géneros Fúngicos (ITS)",   n_top = 30,
                   filename = "heatmap_fungi_genera.pdf")

# ---- 8. Diagrama de Venn: taxa compartidos/únicos entre sitios ----
cat("\n===== DIAGRAMAS DE VENN =====\n")

plot_venn <- function(mat, meta_df, title_str, filename) {
  # Separar muestras por sitio
  site_labels_local <- unique(meta_df[colnames(mat), "site"])
  groups_samples <- split(colnames(mat), meta_df[colnames(mat), "site"])

  # Taxa presentes en cada sitio (suma de muestras del sitio > 0)
  taxa_per_site <- lapply(groups_samples, function(samples) {
    subset_mat <- mat[, samples, drop = FALSE]
    rownames(subset_mat)[rowSums(subset_mat) > 0]
  })

  names(taxa_per_site) <- c("Rancho Gil", "Rancho San Ignacio")

  # Venn con ggVennDiagram
  p <- ggVennDiagram(
    taxa_per_site,
    label_alpha = 0,
    category.names = names(taxa_per_site)
  ) +
    scale_fill_gradient(low = "#e3f2fd", high = "#0d47a1") +
    scale_color_manual(values = c("grey30", "grey30")) +
    labs(title = title_str, fill = "# taxa") +
    theme(
      plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position  = "right"
    )

  ggsave(file.path(OUT_DIR, filename), p, width = 7, height = 5, dpi = 300)
  cat(sprintf("  ✓ Diagrama de Venn guardado: %s\n", filename))

  # Resumen numérico
  shared   <- intersect(taxa_per_site[[1]], taxa_per_site[[2]])
  only_g1  <- setdiff(taxa_per_site[[1]], taxa_per_site[[2]])
  only_g2  <- setdiff(taxa_per_site[[2]], taxa_per_site[[1]])
  cat(sprintf("    Rancho Gil únicos: %d | Compartidos: %d | Rancho San Ignacio únicos: %d\n",
              length(only_g1), length(shared), length(only_g2)))

  return(invisible(list(shared = shared, unique_g1 = only_g1, unique_g2 = only_g2)))
}

venn_bact  <- plot_venn(bact_mat,  meta_clean[colnames(bact_mat), ],
                         "Taxa bacterianos compartidos entre sitios",
                         "venn_bacteria.png")

venn_fungi <- plot_venn(fungi_mat, meta_clean[colnames(fungi_mat), ],
                         "Taxa fúngicos compartidos entre sitios",
                         "venn_fungi.png")

# ---- 9. Bubble chart: Top taxa más abundantes ----
cat("\n===== BUBBLE CHARTS =====\n")

plot_bubble <- function(mat, meta_df, n_top = 15, title_str, filename) {
  top_taxa <- names(sort(rowMeans(mat), decreasing = TRUE))[1:n_top]
  mat_top  <- mat[top_taxa, , drop = FALSE]

  # Formato largo
  df <- as.data.frame(mat_top) %>%
    rownames_to_column("Taxa") %>%
    pivot_longer(-Taxa, names_to = "sample_id", values_to = "abundance") %>%
    left_join(meta_df %>% rownames_to_column("sample_id") %>% select(sample_id, site),
              by = "sample_id") %>%
    mutate(
      Taxa       = factor(Taxa, levels = rev(top_taxa)),
      site_label = site_labels[site]
    ) %>%
    group_by(Taxa, site_label) %>%
    summarise(
      mean_abund = mean(abundance),
      sd_abund   = sd(abundance),
      .groups    = "drop"
    )

  p <- ggplot(df, aes(x = site_label, y = Taxa,
                       size = mean_abund, fill = mean_abund)) +
    geom_point(shape = 21, color = "grey40", alpha = 0.85) +
    scale_size_continuous(range = c(2, 14), name = "Abundancia\nmedia (%)") +
    scale_fill_gradientn(
      colors = c("#e3f2fd", "#1565c0", "#0a47a1"),
      name   = "Abundancia\nmedia (%)"
    ) +
    labs(
      title = title_str,
      x     = "Sitio",
      y     = "Taxa"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y      = element_text(size = 8, face = "italic"),
      axis.text.x      = element_text(size = 10, angle = 15, hjust = 0.8),
      panel.grid.major = element_line(color = "grey92"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11),
      legend.position  = "right"
    )

  ggsave(file.path(OUT_DIR, filename), p, width = 9, height = 8, dpi = 300)
  cat(sprintf("  ✓ Bubble chart guardado: %s\n", filename))

  return(invisible(p))
}

plot_bubble(bact_mat,  meta_clean[colnames(bact_mat), ],
            n_top = 15, title_str = "Top 15 Taxa Bacterianos por Sitio",
            filename = "bubble_bacteria_top15.pdf")

plot_bubble(fungi_mat, meta_clean[colnames(fungi_mat), ],
            n_top = 15, title_str = "Top 15 Taxa Fúngicos por Sitio",
            filename = "bubble_fungi_top15.pdf")

# ---- 10. Tabla: taxa más abundantes con estadísticos descriptivos ----
cat("\n===== EXPORTANDO TABLAS DE COMPOSICIÓN =====\n")

summarise_taxa <- function(mat, meta_df, n_top = 30) {
  groups_col <- meta_df[colnames(mat), "site"]

  site1 <- "Rancho_Gil"
  site2 <- "Rancho_San_Ignacio"

  s1_cols <- colnames(mat)[groups_col == site1]
  s2_cols <- colnames(mat)[groups_col == site2]

  df <- data.frame(
    Taxa             = rownames(mat),
    Mean_overall     = round(rowMeans(mat), 4),
    SD_overall       = round(apply(mat, 1, sd), 4),
    Mean_GIL         = round(rowMeans(mat[, s1_cols, drop = FALSE]), 4),
    SD_GIL           = round(apply(mat[, s1_cols, drop = FALSE], 1, sd), 4),
    Mean_SAN_IGNACIO = round(rowMeans(mat[, s2_cols, drop = FALSE]), 4),
    SD_SAN_IGNACIO   = round(apply(mat[, s2_cols, drop = FALSE], 1, sd), 4),
    Prevalence_total = rowSums(mat > 0),
    Prevalence_GIL   = rowSums(mat[, s1_cols, drop = FALSE] > 0),
    Prevalence_SAN   = rowSums(mat[, s2_cols, drop = FALSE] > 0)
  ) %>%
    arrange(desc(Mean_overall)) %>%
    head(n_top)

  return(df)
}

taxa_summary_bact  <- summarise_taxa(bact_mat,  meta_clean[colnames(bact_mat), ])
taxa_summary_fungi <- summarise_taxa(fungi_mat, meta_clean[colnames(fungi_mat), ])

write.csv(taxa_summary_bact,  file.path(OUT_DIR, "top30_taxa_bacteria.csv"),  row.names = FALSE)
write.csv(taxa_summary_fungi, file.path(OUT_DIR, "top30_taxa_fungi.csv"),     row.names = FALSE)

cat("  ✓ Tablas guardadas: top30_taxa_bacteria.csv, top30_taxa_fungi.csv\n")
cat("\n  Top 10 taxa bacterianos por abundancia media:\n")
print(head(taxa_summary_bact[, c("Taxa", "Mean_overall", "Mean_GIL", "Mean_SAN_IGNACIO")], 10))
cat("\n  Top 10 taxa fúngicos por abundancia media:\n")
print(head(taxa_summary_fungi[, c("Taxa", "Mean_overall", "Mean_GIL", "Mean_SAN_IGNACIO")], 10))

cat("\n✓ Análisis de composición taxonómica completado. Procede con 05_differential_abundance.R\n")

