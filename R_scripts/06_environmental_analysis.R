# =============================================================================
# 06_environmental_analysis.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Análisis de relación entre variables fisicoquímicas del suelo
#              y la composición del microbioma.
#
# Análisis incluidos:
#   1. Correlación entre variables fisicoquímicas (corrplot, heatmap)
#   2. RDA (Redundancy Analysis) - variables → composición microbiana
#   3. CCA (Canonical Correspondence Analysis) - para datos de conteos
#   4. db-RDA (distance-based RDA) con Bray-Curtis
#   5. Mantel test (correlación entre distancias microbianas y ambientales)
#   6. Heatmap de correlación taxa ↔ variables ambientales (Spearman)
#   7. Scatterplots: diversidad alfa vs. variables clave
#
# Dependencias: 01_data_import.R, 02_alpha_diversity.R
# =============================================================================

# ---- 0. Cargar datos y scripts previos ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("01_data_import.R")
source("02_alpha_diversity.R")

suppressPackageStartupMessages({
  library(vegan)
  library(ape)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(corrplot)
  library(RColorBrewer)
  library(patchwork)
  library(pheatmap)
  library(ggrepel)
  library(scales)
})

# ---- 1. Configuración ----
OUT_DIR <- "figures/environmental_analysis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

site_colors <- c("Rancho_Gil" = "#2196F3", "Rancho_San_Ignacio" = "#FF5722")
site_labels <- c("Rancho_Gil" = "Rancho Gil", "Rancho_San_Ignacio" = "Rancho San Ignacio")

# ---- 2. Preparar matriz de variables ambientales ----
# Seleccionar variables numéricas fisicoquímicas
# (excluir columnas de texto y variables derivadas del mapeo)
env_vars <- meta_clean %>%
  select(where(is.numeric)) %>%
  # Eliminar columnas con varianza 0 (constantes entre muestras)
  select(where(~ var(., na.rm = TRUE) > 0)) %>%
  # Eliminar columnas que son exactamente iguales (duplicadas)
  as.matrix()

cat("\n===== VARIABLES AMBIENTALES =====\n")
cat(sprintf("  Variables seleccionadas: %d\n", ncol(env_vars)))
cat("  Nombres:", paste(colnames(env_vars), collapse = ", "), "\n")

# Verificar valores faltantes
na_count <- colSums(is.na(env_vars))
if (any(na_count > 0)) {
  cat("  ⚠ Columnas con NA:\n")
  print(na_count[na_count > 0])
  # Imputar NA con la media de la columna
  for (col in names(na_count[na_count > 0])) {
    env_vars[is.na(env_vars[, col]), col] <- mean(env_vars[, col], na.rm = TRUE)
  }
}

# Estandarizar (z-score) para RDA/CCA
env_scaled <- scale(env_vars)
cat("  ✓ Variables estandarizadas (z-score)\n")

# Asegurar que el orden de muestras coincide
stopifnot(all(rownames(env_vars) == colnames(bact_mat)))

# ---- 3. Correlaciones entre variables ambientales ----
cat("\n===== CORRELACIONES ENTRE VARIABLES AMBIENTALES =====\n")

# Matriz de correlación de Pearson
env_cor <- cor(env_vars, method = "pearson", use = "pairwise.complete.obs")

# p-valores de correlación
env_cor_pval <- matrix(NA, ncol(env_vars), ncol(env_vars))
for (i in 1:ncol(env_vars)) {
  for (j in 1:ncol(env_vars)) {
    ct <- tryCatch(
      cor.test(env_vars[, i], env_vars[, j], method = "pearson"),
      error = function(e) list(p.value = 1)
    )
    env_cor_pval[i, j] <- ct$p.value
  }
}
rownames(env_cor_pval) <- colnames(env_cor_pval) <- colnames(env_vars)

# Corrplot
pdf(file.path(OUT_DIR, "corrplot_environmental_vars.pdf"), width = 10, height = 10)
corrplot(
  env_cor,
  method     = "color",
  type       = "upper",
  order      = "hclust",
  hclust.method = "ward.D2",
  tl.col     = "black",
  tl.srt     = 45,
  tl.cex     = 0.7,
  cl.cex     = 0.8,
  addCoef.col = "black",
  number.cex = 0.5,
  col        = colorRampPalette(c("#2166ac", "white", "#d73027"))(200),
  title      = "Correlación entre variables fisicoquímicas del suelo",
  mar        = c(0, 0, 2, 0),
  p.mat      = env_cor_pval,
  sig.level  = 0.05,
  insig      = "blank"  # ocultar correlaciones no significativas
)
dev.off()
cat("  ✓ Corrplot guardado: corrplot_environmental_vars.pdf\n")

# ---- 4. RDA (Redundancy Analysis) ----
# RDA: la composición microbiana es la variable respuesta;
# las variables ambientales son los predictores.
cat("\n===== RDA (Redundancy Analysis) =====\n")

# Transformación Hellinger de la matriz de abundancia
bact_hell  <- decostand(t(bact_mat),  method = "hellinger")
fungi_hell <- decostand(t(fungi_mat), method = "hellinger")

# RDA con todas las variables (puede estar sobre-ajustada con n=6)
rda_bact  <- rda(bact_hell  ~ ., data = as.data.frame(env_scaled))
rda_fungi <- rda(fungi_hell ~ ., data = as.data.frame(env_scaled))

cat("\n  RDA Bacterias:\n"); print(summary(rda_bact, display = NULL))
cat("\n  RDA Hongos:\n");    print(summary(rda_fungi, display = NULL))

# R² ajustado (medida de ajuste del modelo)
r2_bact_adj  <- RsquareAdj(rda_bact)$adj.r.squared
r2_fungi_adj <- RsquareAdj(rda_fungi)$adj.r.squared
cat(sprintf("\n  R² ajustado RDA Bacterias: %.4f\n", r2_bact_adj))
cat(sprintf("  R² ajustado RDA Hongos:    %.4f\n", r2_fungi_adj))

# Test de significancia global (permutaciones)
set.seed(42)
perm_rda_bact  <- anova.cca(rda_bact,  permutations = 999, by = "axis")
perm_rda_fungi <- anova.cca(rda_fungi, permutations = 999, by = "axis")
cat("\n  ANOVA-RDA ejes - Bacterias:\n"); print(perm_rda_bact)
cat("\n  ANOVA-RDA ejes - Hongos:\n");    print(perm_rda_fungi)

# Selección de variables por pasos (forward selection)
# Útil para identificar las variables más importantes con n pequeño
cat("\n  Forward selection de variables (bacterias):\n")
rda_null_bact <- rda(bact_hell ~ 1, data = as.data.frame(env_scaled))
set.seed(42)
rda_fwd_bact  <- tryCatch(
  ordiR2step(rda_null_bact, scope = formula(rda_bact),
             direction = "forward", R2permutations = 999, R2scope = TRUE),
  error = function(e) {
    cat("  ⚠ Forward selection fallido (n muy pequeño):", e$message, "\n")
    NULL
  }
)
if (!is.null(rda_fwd_bact)) print(rda_fwd_bact)

# ---- 5. Función: Triplot RDA/CCA ----
plot_rda_triplot <- function(rda_obj, meta_df, groups_vec, title_str, filename) {
  sc <- scores(rda_obj, display = c("sites", "bp"), scaling = 2)

  # Coordenadas de sitios
  sites_df <- as.data.frame(sc$sites)
  names(sites_df)[1:2] <- c("RDA1", "RDA2")
  sites_df$sample_id   <- rownames(sites_df)
  sites_df$site        <- groups_vec

  # Coordenadas de biplot (flechas de variables ambientales)
  bp_df <- as.data.frame(sc$biplot)
  names(bp_df)[1:2] <- c("RDA1", "RDA2")
  bp_df$Variable <- rownames(bp_df)

  # Escala de flechas (ajustar al espacio del gráfico)
  arrow_scale <- 0.5 * max(abs(range(sites_df[, 1:2])))
  bp_df[, 1:2] <- bp_df[, 1:2] * arrow_scale

  # Porcentaje de varianza explicado por cada eje
  var_exp <- eigenvals(rda_obj) / sum(eigenvals(rda_obj))
  xlab <- sprintf("RDA1 (%.1f%%)", var_exp[1] * 100)
  ylab <- sprintf("RDA2 (%.1f%%)", var_exp[2] * 100)

  p <- ggplot() +
    # Flechas de variables ambientales
    geom_segment(data = bp_df,
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 color = "grey40", linewidth = 0.5) +
    geom_text_repel(data = bp_df,
                    aes(x = RDA1, y = RDA2, label = Variable),
                    size = 2.5, color = "grey25", fontface = "italic",
                    max.overlaps = 20) +
    # Puntos de sitios
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
    geom_point(data = sites_df,
               aes(x = RDA1, y = RDA2, color = site, shape = site),
               size = 5, alpha = 0.9) +
    geom_text_repel(data = sites_df,
                    aes(x = RDA1, y = RDA2, label = sample_id, color = site),
                    size = 3, nudge_y = 0.02, show.legend = FALSE) +
    scale_color_manual(values = site_colors, labels = site_labels, name = "Sitio") +
    scale_shape_manual(values = c(16, 17), labels = site_labels, name = "Sitio") +
    labs(title = title_str,
         subtitle = sprintf("R² ajustado = %.3f", RsquareAdj(rda_obj)$adj.r.squared),
         x = xlab, y = ylab) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11)
    )

  ggsave(file.path(OUT_DIR, filename), p, width = 9, height = 7, dpi = 300)
  cat(sprintf("  ✓ Triplot RDA guardado: %s\n", filename))
  return(invisible(p))
}

groups <- meta_clean[colnames(bact_mat), "site"]
plot_rda_triplot(rda_bact,  meta_clean, groups,
                 "RDA Bacterias/Arqueas vs. Variables Ambientales",
                 "rda_triplot_bacteria.pdf")
plot_rda_triplot(rda_fungi, meta_clean, groups,
                 "RDA Hongos vs. Variables Ambientales",
                 "rda_triplot_fungi.pdf")

# ---- 6. db-RDA (distance-based RDA) ----
cat("\n===== db-RDA (Bray-Curtis) =====\n")

dist_bray_bact  <- vegdist(t(bact_mat),  method = "bray")
dist_bray_fungi <- vegdist(t(fungi_mat), method = "bray")

dbrda_bact  <- dbrda(dist_bray_bact  ~ ., data = as.data.frame(env_scaled))
dbrda_fungi <- dbrda(dist_bray_fungi ~ ., data = as.data.frame(env_scaled))

cat("  db-RDA Bacterias R²:\n"); print(RsquareAdj(dbrda_bact))
cat("  db-RDA Hongos R²:\n");    print(RsquareAdj(dbrda_fungi))

# ---- 7. Mantel Test ----
cat("\n===== MANTEL TEST =====\n")
cat("  Evalúa si la similitud ambiental predice la similitud microbiana.\n\n")

# Distancia ambiental (Euclidean sobre variables estandarizadas)
dist_env <- dist(env_scaled, method = "euclidean")

# Mantel test
set.seed(42)
mantel_bact  <- mantel(dist_bray_bact,  dist_env, method = "spearman",
                        permutations = 999, na.rm = TRUE)
mantel_fungi <- mantel(dist_bray_fungi, dist_env, method = "spearman",
                        permutations = 999, na.rm = TRUE)

cat(sprintf("  Mantel Bacterias: r = %.4f | p = %.4f\n",
            mantel_bact$statistic, mantel_bact$signif))
cat(sprintf("  Mantel Hongos:    r = %.4f | p = %.4f\n",
            mantel_fungi$statistic, mantel_fungi$signif))

# Mantel parcial (controlando por variable de grupo)
# Útil para separar el efecto del sitio del efecto de las variables continuas
site_dist <- vegdist(model.matrix(~ groups - 1), method = "euclidean")
set.seed(42)
mantel_part_bact <- mantel.partial(dist_bray_bact,  dist_env, site_dist,
                                    permutations = 999)
mantel_part_fungi <- mantel.partial(dist_bray_fungi, dist_env, site_dist,
                                     permutations = 999)
cat(sprintf("  Mantel parcial Bacterias (controlando sitio): r = %.4f | p = %.4f\n",
            mantel_part_bact$statistic, mantel_part_bact$signif))
cat(sprintf("  Mantel parcial Hongos:                        r = %.4f | p = %.4f\n",
            mantel_part_fungi$statistic, mantel_part_fungi$signif))

# ---- 8. Correlación taxa ↔ variables ambientales (Spearman) ----
cat("\n===== CORRELACIÓN TAXA vs. VARIABLES AMBIENTALES =====\n")

compute_taxa_env_cor <- function(mat, env_mat, n_top = 30) {
  # Seleccionar top taxa por abundancia media
  top_taxa <- names(sort(rowMeans(mat), decreasing = TRUE))[1:min(n_top, nrow(mat))]
  mat_top  <- mat[top_taxa, , drop = FALSE]

  # Calcular correlación de Spearman entre cada taxón y cada variable ambiental
  cor_mat  <- matrix(NA, nrow = length(top_taxa), ncol = ncol(env_mat))
  pval_mat <- matrix(NA, nrow = length(top_taxa), ncol = ncol(env_mat))
  rownames(cor_mat)  <- rownames(pval_mat) <- top_taxa
  colnames(cor_mat)  <- colnames(pval_mat) <- colnames(env_mat)

  for (taxon in top_taxa) {
    for (var in colnames(env_mat)) {
      ct <- tryCatch(
        cor.test(as.numeric(mat_top[taxon, ]), env_mat[, var],
                 method = "spearman", exact = FALSE),
        error = function(e) list(estimate = NA, p.value = NA)
      )
      cor_mat[taxon, var]  <- ct$estimate
      pval_mat[taxon, var] <- ct$p.value
    }
  }

  return(list(cor = cor_mat, pval = pval_mat))
}

taxa_env_bact  <- compute_taxa_env_cor(bact_mat,  env_vars, n_top = 30)
taxa_env_fungi <- compute_taxa_env_cor(fungi_mat, env_vars, n_top = 30)

# Heatmap de correlación taxa-ambiente
plot_taxa_env_heatmap <- function(cor_list, title_str, filename) {
  cor_mat  <- cor_list$cor
  pval_mat <- cor_list$pval

  # Crear matriz de asteriscos de significancia
  sig_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
  sig_mat[pval_mat < 0.05]  <- "*"
  sig_mat[pval_mat < 0.01]  <- "**"
  sig_mat[pval_mat < 0.001] <- "***"
  sig_mat[is.na(pval_mat)]  <- ""

  cor_mat[is.na(cor_mat)] <- 0

  pdf(file.path(OUT_DIR, filename), width = 12, height = 10)
  pheatmap(
    cor_mat,
    color             = colorRampPalette(c("#2166ac", "white", "#d73027"))(100),
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    clustering_method = "ward.D2",
    display_numbers   = sig_mat,
    number_color      = "black",
    fontsize_number   = 9,
    fontsize_row      = 7,
    fontsize_col      = 8,
    border_color      = "grey90",
    main              = paste(title_str, "\n(* p<0.05, ** p<0.01, *** p<0.001)"),
    angle_col         = 45
  )
  dev.off()
  cat(sprintf("  ✓ Heatmap correlación taxa-ambiente guardado: %s\n", filename))
}

plot_taxa_env_heatmap(taxa_env_bact,  "Correlación Spearman: Top 30 Bacterias × Variables Ambientales",
                       "heatmap_taxa_env_bacteria.pdf")
plot_taxa_env_heatmap(taxa_env_fungi, "Correlación Spearman: Top 30 Hongos × Variables Ambientales",
                       "heatmap_taxa_env_fungi.pdf")

# ---- 9. Scatterplots: diversidad alfa vs. variables ambientales ----
cat("\n===== SCATTERPLOTS: DIVERSIDAD ALFA vs. VARIABLES AMBIENTALES =====\n")

# Unir diversidad alfa con variables ambientales
alpha_env_bact <- alpha_bact %>%
  column_to_rownames("sample_id") %>%
  bind_cols(as.data.frame(env_vars)[rownames(.), ])

# Variables más relevantes para mostrar (selección manual)
key_env_vars <- c("ph", "materia_organica_", "humedad_", "conductividad_elctrica_ms_cm",
                   "nitrгeno_inorgnico_mgkg")

# Ajustar nombres según los que existan en los datos
key_env_vars <- key_env_vars[key_env_vars %in% colnames(alpha_env_bact)]

if (length(key_env_vars) > 0) {
  plot_list <- lapply(key_env_vars, function(var) {
    ggplot(alpha_env_bact, aes(x = .data[[var]], y = Shannon, color = site)) +
      geom_point(size = 4, alpha = 0.9) +
      geom_smooth(method = "lm", se = TRUE, color = "grey40",
                  fill = "grey80", linewidth = 0.5) +
      scale_color_manual(values = site_colors, labels = site_labels) +
      ggpubr::stat_cor(method = "spearman", size = 3) +
      labs(x = var, y = "Shannon H'", color = "Sitio") +
      theme_bw(base_size = 10) +
      theme(legend.position = "bottom")
  })

  fig_scatter <- wrap_plots(plot_list, ncol = min(3, length(plot_list))) +
    plot_annotation(
      title = "Diversidad alfa (Shannon) vs. Variables Ambientales — Bacterias",
      theme = theme(plot.title = element_text(face = "bold"))
    )

  ggsave(file.path(OUT_DIR, "scatterplots_alpha_vs_env.pdf"),
         fig_scatter, width = 12, height = 5, dpi = 300)
  cat("  ✓ Scatterplots guardados: scatterplots_alpha_vs_env.pdf\n")
}

# ---- 10. Exportar resultados ----
cat("\n===== EXPORTANDO RESULTADOS =====\n")

mantel_summary <- data.frame(
  Dataset   = c("Bacterias", "Hongos", "Bacterias (parcial)", "Hongos (parcial)"),
  r_Mantel  = c(mantel_bact$statistic, mantel_fungi$statistic,
                 mantel_part_bact$statistic, mantel_part_fungi$statistic),
  p_Mantel  = c(mantel_bact$signif, mantel_fungi$signif,
                 mantel_part_bact$signif, mantel_part_fungi$signif)
)
write.csv(mantel_summary, file.path(OUT_DIR, "mantel_test_results.csv"), row.names = FALSE)

write.csv(as.data.frame(taxa_env_bact$cor),
          file.path(OUT_DIR, "spearman_cor_taxa_env_bacteria.csv"))
write.csv(as.data.frame(taxa_env_fungi$cor),
          file.path(OUT_DIR, "spearman_cor_taxa_env_fungi.csv"))

cat("  ✓ Tablas exportadas: mantel_test_results.csv, spearman_cor_taxa_env_*.csv\n")

cat("\n✓ Análisis de variables ambientales completado. Procede con 07_network_analysis.R\n")
