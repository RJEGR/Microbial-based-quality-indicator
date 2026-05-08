# =============================================================================
# 03_beta_diversity.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Análisis de diversidad beta (diferencias entre comunidades).
#
# Análisis incluidos:
#   1. Distancias de disimilitud (Bray-Curtis, Jaccard, Euclidean)
#   2. NMDS (Non-metric Multidimensional Scaling)
#   3. PCoA (Principal Coordinates Analysis)
#   4. PCA (Principal Component Analysis sobre datos transformados)
#   5. Clustering jerárquico (hclust)
#   6. PERMANOVA (adonis2) - prueba multivariante de diferencias
#   7. ANOSIM - análisis de similitud
#   8. Análisis de homogeneidad de dispersión (betadisper)
#
# Dependencias: 01_data_import.R
# =============================================================================

# ---- 0. Cargar datos ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("01_data_import.R")

suppressPackageStartupMessages({
  library(vegan)
  library(ape)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
})

# ---- 1. Configuración ----

OUT_DIR1 <- here::here("~/Documents/Metagenomics/Microbial-soil/Artículo/")  # Ajustar si los scripts están en R_scripts/

OUT_DIR2 <- "figures/beta_diversity"

OUT_DIR <- file.path(OUT_DIR1, OUT_DIR2)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

site_colors <- c("Rancho_Gil" = "#2196F3", "Rancho_San_Ignacio" = "#FF5722")
site_shapes <- c("Rancho_Gil" = 16,        "Rancho_San_Ignacio" = 17)
site_labels <- c("Rancho_Gil" = "Rancho Gil", "Rancho_San_Ignacio" = "Rancho San Ignacio")

# Grupos de muestras
groups <- meta_clean[colnames(bact_mat), "site"]

# ---- 2. Transformación de datos ----
# Para análisis de composición, se recomienda transformación de Hellinger o
# transformación de raíz cuadrada antes de PCA/RDA
bact_hell  <- decostand(t(bact_mat),  method = "hellinger")  # Hellinger: buen default
fungi_hell <- decostand(t(fungi_mat), method = "hellinger")

# CLR (Centered Log-Ratio) para análisis composicional robusto
# Requiere añadir pseudocount pequeño para evitar log(0)
clr_transform <- function(mat_t) {
  mat_t[mat_t == 0] <- 1e-6  # pseudocount mínimo
  clr <- log(mat_t) - rowMeans(log(mat_t))
  return(clr)
}
bact_clr  <- clr_transform(t(bact_mat))
fungi_clr <- clr_transform(t(fungi_mat))

# ---- 3. Matrices de distancia ----
# Bray-Curtis (estándar en microbioma: semi-métrica, buena para datos composicionales)
dist_bray_bact  <- vegdist(t(bact_mat),  method = "bray")
dist_bray_fungi <- vegdist(t(fungi_mat), method = "bray")

# Jaccard (basada en presencia/ausencia)
dist_jac_bact  <- vegdist(t(bact_mat),  method = "jaccard", binary = TRUE)
dist_jac_fungi <- vegdist(t(fungi_mat), method = "jaccard", binary = TRUE)

# Euclidean sobre Hellinger (equivalente a distancia de Chord)
dist_hell_bact  <- dist(bact_hell)
dist_hell_fungi <- dist(fungi_hell)

cat("  ✓ Matrices de distancia calculadas\n")
cat("\n  Bray-Curtis - Bacterias:\n"); print(round(as.matrix(dist_bray_bact), 3))
cat("\n  Bray-Curtis - Hongos:\n");    print(round(as.matrix(dist_bray_fungi), 3))

# ---- 4. PERMANOVA (adonis2) ----
cat("\n===== PERMANOVA (adonis2) =====\n")
cat("  ADVERTENCIA: Con n=6 total (3 por grupo), el número de permutaciones\n")
cat("  posibles es muy limitado (~10 permutaciones exactas). Los p-valores\n")
cat("  mínimos obtenibles son p≈0.05 (1 de 10 permutaciones posibles).\n\n")

# Función para ejecutar PERMANOVA con reporte
run_permanova <- function(dist_mat, meta_df, formula_str = "site",
                          nperm = 999, label = "") {
  set.seed(42)  # Reproducibilidad
  result <- adonis2(
    as.formula(paste("dist_mat ~", formula_str)),
    data = meta_df,
    permutations = nperm,
    method = "bray"  # No aplica cuando ya se provee dist_mat, pero requerido
  )
  cat(sprintf("\n  PERMANOVA - %s:\n", label))
  print(result)
  return(result)
}

perm_bact_bray  <- run_permanova(dist_bray_bact,  meta_clean[colnames(bact_mat), ],
                                  label = "Bacterias (Bray-Curtis)")
perm_fungi_bray <- run_permanova(dist_bray_fungi, meta_clean[colnames(fungi_mat), ],
                                  label = "Hongos (Bray-Curtis)")

# ---- 5. ANOSIM ----
cat("\n===== ANOSIM =====\n")
set.seed(42)
anosim_bact  <- anosim(dist_bray_bact,  grouping = groups, permutations = 999)
anosim_fungi <- anosim(dist_bray_fungi, grouping = groups, permutations = 999)

cat("\n  ANOSIM - Bacterias:\n")
cat(sprintf("    R = %.4f | p = %.4f\n", anosim_bact$statistic, anosim_bact$signif))

cat("\n  ANOSIM - Hongos:\n")
cat(sprintf("    R = %.4f | p = %.4f\n", anosim_fungi$statistic, anosim_fungi$signif))

# ---- 6. Homogeneidad de dispersión (betadisper) ----
cat("\n===== BETADISPER (homogeneidad de varianzas) =====\n")
betad_bact  <- betadisper(dist_bray_bact,  group = groups)
betad_fungi <- betadisper(dist_bray_fungi, group = groups)

cat("\n  Dispersión media por grupo (Bacterias):\n")
print(betad_bact$group.distances)
cat("\n  Dispersión media por grupo (Hongos):\n")
print(betad_fungi$group.distances)

# Test de Levene generalizado (permutest)
set.seed(42)
perm_disp_bact  <- permutest(betad_bact,  permutations = 99)
perm_disp_fungi <- permutest(betad_fungi, permutations = 99)
cat("\n  p-valor betadisper Bacterias:", round(perm_disp_bact$tab["Groups","Pr(>F)"], 4), "\n")
cat("  p-valor betadisper Hongos:   ", round(perm_disp_fungi$tab["Groups","Pr(>F)"], 4), "\n")

# ---- 7. Función genérica: extraer coordenadas de ordenación ----
extract_ordination <- function(ord_obj, type = c("nmds", "pcoa", "pca"),
                                meta_df = NULL, groups_vec = NULL) {
  type <- match.arg(type)

  if (type == "nmds") {
    scores <- as.data.frame(scores(ord_obj, display = "sites"))
    names(scores)[1:2] <- c("Axis1", "Axis2")
    stress <- ord_obj$stress
    xlab <- paste0("NMDS1")
    ylab <- paste0("NMDS2")
    subtitle <- sprintf("Stress = %.4f (< 0.2 = bueno)", stress)
  } else if (type == "pcoa") {
    scores <- as.data.frame(ord_obj$vectors[, 1:2])
    names(scores)[1:2] <- c("Axis1", "Axis2")
    pct <- round(ord_obj$values$Relative_eig[1:2] * 100, 1)
    xlab <- paste0("PCoA1 (", pct[1], "%)")
    ylab <- paste0("PCoA2 (", pct[2], "%)")
    subtitle <- ""
  } else if (type == "pca") {
    scores <- as.data.frame(ord_obj$x[, 1:2])
    names(scores)[1:2] <- c("Axis1", "Axis2")
    var_exp <- summary(ord_obj)$importance
    pct <- round(var_exp["Proportion of Variance", 1:2] * 100, 1)
    xlab <- paste0("PC1 (", pct[1], "%)")
    ylab <- paste0("PC2 (", pct[2], "%)")
    subtitle <- ""
  }

  scores$sample_id <- rownames(scores)
  if (!is.null(groups_vec)) scores$site <- groups_vec
  if (!is.null(meta_df)) {
    scores <- scores %>%
      left_join(meta_df %>% rownames_to_column("sample_id"), by = "sample_id")
  }

  attr(scores, "xlab")     <- xlab
  attr(scores, "ylab")     <- ylab
  attr(scores, "subtitle") <- subtitle

  return(scores)
}

# ---- 8. NMDS ----
cat("\n===== NMDS (Non-metric Multidimensional Scaling) =====\n")
set.seed(42)
nmds_bact  <- metaMDS(dist_bray_bact,  k = 2, trymax = 100, trace = FALSE)
nmds_fungi <- metaMDS(dist_bray_fungi, k = 2, trymax = 100, trace = FALSE)

cat(sprintf("  Stress NMDS Bacterias: %.4f\n", nmds_bact$stress))
cat(sprintf("  Stress NMDS Hongos:    %.4f\n", nmds_fungi$stress))
if (nmds_bact$stress > 0.2)  cat("  ⚠ Stress alto en bacterias (> 0.2): ordination poco confiable.\n")
if (nmds_fungi$stress > 0.2) cat("  ⚠ Stress alto en hongos (> 0.2): ordination poco confiable.\n")

nmds_bact_df  <- extract_ordination(nmds_bact,  type = "nmds", groups_vec = groups)
nmds_fungi_df <- extract_ordination(nmds_fungi, type = "nmds", groups_vec = groups)

# ---- 9. PCoA ----
cat("\n===== PCoA (Principal Coordinates Analysis) =====\n")
pcoa_bact  <- pcoa(dist_bray_bact)
pcoa_fungi <- pcoa(dist_bray_fungi)

pcoa_bact_df  <- extract_ordination(pcoa_bact,  type = "pcoa", groups_vec = groups)
pcoa_fungi_df <- extract_ordination(pcoa_fungi, type = "pcoa", groups_vec = groups)

# ---- 10. PCA sobre datos transformados (Hellinger) ----
cat("\n===== PCA (sobre transformación Hellinger) =====\n")
pca_bact  <- prcomp(bact_hell,  center = TRUE, scale. = FALSE)
pca_fungi <- prcomp(fungi_hell, center = TRUE, scale. = FALSE)

pca_bact_df  <- extract_ordination(pca_bact,  type = "pca", groups_vec = groups)
pca_fungi_df <- extract_ordination(pca_fungi, type = "pca", groups_vec = groups)

# ---- 11. Función de visualización de ordenación ----
plot_ordination <- function(df, x = "Axis1", y = "Axis2",
                             xlab = "Eje 1", ylab = "Eje 2",
                             title = "", subtitle = "",
                             label_col = "sample_id") {
  # Extraer atributos si están disponibles
  if (!is.null(attr(df, "xlab")))     xlab     <- attr(df, "xlab")
  if (!is.null(attr(df, "ylab")))     ylab     <- attr(df, "ylab")
  if (!is.null(attr(df, "subtitle"))) subtitle <- attr(df, "subtitle")

  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]],
                       color = site, shape = site, label = .data[[label_col]])) +
    # Elipses de confianza 95% (sólo informativas con n=3)
    stat_ellipse(aes(fill = site), geom = "polygon", alpha = 0.08,
                 level = 0.95, type = "t", linetype = 2, linewidth = 0.4) +
    geom_point(size = 5, alpha = 0.9) +
    geom_text(nudge_y = 0.02, size = 3, color = "grey30", fontface = "italic") +
    scale_color_manual(values = site_colors, labels = site_labels) +
    scale_fill_manual(values  = site_colors, labels = site_labels) +
    scale_shape_manual(values = site_shapes, labels = site_labels) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    labs(title    = title,
         subtitle = subtitle,
         x        = xlab,
         y        = ylab,
         color    = "Sitio",
         fill     = "Sitio",
         shape    = "Sitio") +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 9, color = "grey40")
    )

  return(p)
}

# ---- 12. Generar y guardar figuras ----
cat("\n===== GENERANDO FIGURAS DE ORDENACIÓN =====\n")

# NMDS
p_nmds_bact  <- plot_ordination(nmds_bact_df,  title = "NMDS - Bacterias/Arqueas (Bray-Curtis)")
p_nmds_fungi <- plot_ordination(nmds_fungi_df, title = "NMDS - Hongos (Bray-Curtis)")

# PCoA
p_pcoa_bact  <- plot_ordination(pcoa_bact_df,  title = "PCoA - Bacterias/Arqueas (Bray-Curtis)")
p_pcoa_fungi <- plot_ordination(pcoa_fungi_df, title = "PCoA - Hongos (Bray-Curtis)")

# PCA
p_pca_bact  <- plot_ordination(pca_bact_df,   title = "PCA - Bacterias/Arqueas (Hellinger)")
p_pca_fungi <- plot_ordination(pca_fungi_df,  title = "PCA - Hongos (Hellinger)")

# Panel 2x3: todos los métodos
fig_ordination <- (p_nmds_bact | p_pcoa_bact | p_pca_bact) /
                  (p_nmds_fungi | p_pcoa_fungi | p_pca_fungi) +
  plot_annotation(
    title    = "Diversidad Beta — Suelo de Vid (Vitis vinifera)",
    subtitle = "Fila superior: Bacterias/Arqueas (16S) | Fila inferior: Hongos (ITS)",
    caption  = "Elipses de confianza 95% (informativas, n=3). PERMANOVA ver estadísticos.",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      plot.caption  = element_text(size = 8,  color = "grey50")
    )
  )

ggsave(file.path(OUT_DIR, "beta_diversity_ordination.pdf"),
       fig_ordination, width = 16, height = 10, dpi = 300)
ggsave(file.path(OUT_DIR, "beta_diversity_ordination.png"),
       fig_ordination, width = 16, height = 10, dpi = 300)
cat("  ✓ Figura guardada: beta_diversity_ordination.pdf/png\n")

# ---- 13. Clustering jerárquico ----
cat("\n===== CLUSTERING JERÁRQUICO =====\n")

plot_hclust <- function(dist_mat, groups_vec, title_str, col_vec) {
  hc <- hclust(dist_mat, method = "ward.D2")

  # Colorear etiquetas por sitio
  label_colors <- col_vec[groups_vec[hc$order]]

  # Usar dendextend para un dendrograma más elegante
  if (requireNamespace("dendextend", quietly = TRUE)) {
    library(dendextend)
    dend <- as.dendrogram(hc)
    dend <- color_labels(dend, col = label_colors)
    dend <- set(dend, "labels_cex", 0.9)

    pdf(file.path(OUT_DIR, paste0("hclust_", gsub(" |/", "_", title_str), ".pdf")),
        width = 8, height = 5)
    par(mar = c(4, 4, 3, 1))
    plot(dend,
         main = title_str,
         sub  = "Método: Ward D2 | Distancia: Bray-Curtis",
         ylab = "Altura (disimilitud)")
    legend("topright", legend = names(col_vec), fill = col_vec,
           bty = "n", cex = 0.9, title = "Sitio")
    dev.off()
  } else {
    pdf(file.path(OUT_DIR, paste0("hclust_", gsub(" |/", "_", title_str), ".pdf")),
        width = 8, height = 5)
    par(mar = c(4, 4, 3, 1))
    plot(hc, main = title_str,
         sub  = "Método: Ward D2 | Distancia: Bray-Curtis",
         ylab = "Altura", xlab = "Muestras")
    dev.off()
  }
  cat(sprintf("  ✓ Dendrograma guardado: hclust_%s.pdf\n", title_str))
}

plot_hclust(dist_bray_bact,  groups, "Bacterias Bray-Curtis", site_colors)
plot_hclust(dist_bray_fungi, groups, "Hongos Bray-Curtis",    site_colors)

# ---- 14. Heatmap de distancias ----
cat("\n===== HEATMAP DE DISTANCIAS =====\n")

plot_dist_heatmap <- function(dist_mat, meta_df, title_str) {
  mat <- as.matrix(dist_mat)
  annot <- data.frame(Sitio = meta_df[rownames(mat), "site"])
  rownames(annot) <- rownames(mat)

  ann_colors <- list(Sitio = site_colors)

  pdf(file.path(OUT_DIR, paste0("heatmap_dist_", gsub(" |/", "_", title_str), ".pdf")),
      width = 6, height = 5)
  pheatmap::pheatmap(
    mat,
    color            = colorRampPalette(c("#1a237e", "#e3f2fd", "#b71c1c"))(50),
    annotation_row   = annot,
    annotation_col   = annot,
    annotation_colors = ann_colors,
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    main             = paste("Distancias Bray-Curtis —", title_str),
    fontsize         = 10,
    display_numbers  = TRUE,
    number_format    = "%.2f",
    number_color     = "white"
  )
  dev.off()
  cat(sprintf("  ✓ Heatmap de distancias guardado: %s\n", title_str))
}

plot_dist_heatmap(dist_bray_bact,  meta_clean[colnames(bact_mat), ],  "Bacterias")
plot_dist_heatmap(dist_bray_fungi, meta_clean[colnames(fungi_mat), ], "Hongos")

# ---- 15. Resumen estadístico ----
cat("\n===== RESUMEN ESTADÍSTICO DIVERSIDAD BETA =====\n")

stats_beta <- data.frame(
  Dataset    = c("Bacterias (Bray-Curtis)", "Hongos (Bray-Curtis)"),
  PERMANOVA_R2 = c(
    round(perm_bact_bray$R2[1], 4),
    round(perm_fungi_bray$R2[1], 4)
  ),
  PERMANOVA_p = c(
    round(perm_bact_bray$`Pr(>F)`[1], 4),
    round(perm_fungi_bray$`Pr(>F)`[1], 4)
  ),
  ANOSIM_R = c(
    round(anosim_bact$statistic, 4),
    round(anosim_fungi$statistic, 4)
  ),
  ANOSIM_p = c(
    round(anosim_bact$signif, 4),
    round(anosim_fungi$signif, 4)
  ),
  NMDS_stress = c(
    round(nmds_bact$stress, 4),
    round(nmds_fungi$stress, 4)
  )
)

print(stats_beta)
write.csv(stats_beta, file.path(OUT_DIR, "beta_diversity_statistics.csv"), row.names = FALSE)
cat("  ✓ Tabla exportada: beta_diversity_statistics.csv\n")

cat("\n✓ Análisis de diversidad beta completado. Procede con 04_taxonomic_composition.R\n")

