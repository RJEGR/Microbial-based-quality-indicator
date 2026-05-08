# =============================================================================
# 02_alpha_diversity.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Análisis de diversidad alfa (diversidad dentro de cada muestra).
#
# Análisis incluidos:
#   1. Cálculo de métricas de diversidad alfa:
#      - Richness (riqueza observada de especies/taxa)
#      - Shannon (H') - diversidad que pondera abundancia
#      - Simpson (1-D) - dominancia/uniformidad
#      - Chao1 - riqueza estimada (para pseudo-conteos)
#      - ACE  - riqueza estimada con corrección de taxa raros
#      - Pielou's Evenness (J) - equitatividad
#   2. Curvas de rarefacción
#   3. Visualización con boxplots (ggplot2 + ggpubr)
#   4. Pruebas estadísticas (Kruskal-Wallis / Wilcoxon Mann-Whitney)
#      NOTA: Con n=3 por grupo, el poder estadístico es MUY limitado.
#
# Dependencias: 01_data_import.R
# =============================================================================

rm(list = ls())

if(!is.null(dev.list())) dev.off()

# ---- 0. Cargar datos (ejecuta el script de importación) ----
# Ajustar la ruta según la ubicación del script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("01_data_import.R")

# Paquetes adicionales para este script
suppressPackageStartupMessages({
  library(vegan)
  library(phyloseq)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  library(patchwork)
  library(scales)
})

# ---- 1. Configuración de salida ----

OUT_DIR <- here::here("~/Documents/Metagenomics/Microbial-soil/Artículo/figures/alpha_diversity/")  # Ajustar si los scripts están en R_scripts/


dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Paleta de colores para los dos sitios
site_colors <- c("Rancho_Gil" = "#2196F3", "Rancho_San_Ignacio" = "#FF5722")
site_labels <- c("Rancho_Gil" = "Rancho Gil", "Rancho_San_Ignacio" = "Rancho San Ignacio")

# ---- 2. Función: calcular diversidad alfa ----
calc_alpha_diversity <- function(mat, pseudo_mat = NULL, group_var = NULL) {
  # mat: matriz de abundancias relativas (taxa en filas, muestras en columnas)
  # pseudo_mat: pseudo-conteos (para Chao1 y ACE)
  # group_var: vector de grupos para cada muestra

  mat_t <- t(mat)  # vegan espera muestras en filas, taxa en columnas

  # Métricas básicas (vegan)
  richness  <- rowSums(mat_t > 0)                     # Riqueza observada
  shannon   <- diversity(mat_t, index = "shannon")    # H' (log natural)
  simpson   <- diversity(mat_t, index = "simpson")    # 1 - D
  evenness  <- shannon / log(richness)                 # Pielou's J (evenness)

  # Métricas de riqueza estimada (requieren conteos enteros → usar pseudo_mat)
  if (!is.null(pseudo_mat)) {
    pseudo_t <- t(pseudo_mat)
    chao1_vals <- estimateR(pseudo_t)  # Retorna: S.obs, S.chao1, se.chao1, S.ACE, se.ACE
    chao1 <- chao1_vals["S.chao1", ]
    ace   <- chao1_vals["S.ACE", ]
  } else {
    chao1 <- rep(NA, ncol(mat))
    ace   <- rep(NA, ncol(mat))
  }

  # Ensamblar data frame
  alpha_df <- data.frame(
    sample_id = colnames(mat),
    Richness  = richness,
    Shannon   = round(shannon, 4),
    Simpson   = round(simpson, 4),
    Evenness  = round(evenness, 4),
    Chao1     = round(chao1, 1),
    ACE       = round(ace, 1),
    stringsAsFactors = FALSE
  )

  # Añadir grupo si se proporciona
  if (!is.null(group_var)) {
    alpha_df$site <- group_var
  }

  return(alpha_df)
}

# ---- 3. Calcular diversidad alfa para bacterias y hongos ----
cat("\n===== DIVERSIDAD ALFA - BACTERIAS/ARQUEAS (16S) =====\n")
alpha_bact <- calc_alpha_diversity(
  mat       = bact_mat,
  pseudo_mat = bact_pseudo,
  group_var  = meta_clean[colnames(bact_mat), "site"]
)
print(alpha_bact)

cat("\n===== DIVERSIDAD ALFA - HONGOS (ITS) =====\n")
alpha_fungi <- calc_alpha_diversity(
  mat       = fungi_mat,
  pseudo_mat = fungi_pseudo,
  group_var  = meta_clean[colnames(fungi_mat), "site"]
)
print(alpha_fungi)

# Añadir etiqueta de reino
alpha_bact$kingdom  <- "Bacteria/Archaea (16S)"
alpha_fungi$kingdom <- "Fungi (ITS)"

# Combinar en un único data frame para visualización
alpha_all <- bind_rows(alpha_bact, alpha_fungi)

# ---- 4. Pruebas estadísticas (Wilcoxon / Kruskal-Wallis) ----
cat("\n===== PRUEBAS ESTADÍSTICAS =====\n")
cat("  ADVERTENCIA: n=3 por grupo. Los p-valores deben interpretarse\n")
cat("  con extrema precaución (muy bajo poder estadístico).\n\n")

metrics <- c("Richness", "Shannon", "Simpson", "Evenness", "Chao1")
kingdoms <- c("Bacteria/Archaea (16S)", "Fungi (ITS)")

stats_results <- list()
for (k in kingdoms) {
  df_k <- alpha_all %>% filter(kingdom == k)
  for (metric in metrics) {
    # Wilcoxon Mann-Whitney (no paramétrico, apropiado para n pequeño)
    test_result <- tryCatch(
      wilcox.test(
        x       = df_k[[metric]][df_k$site == "Rancho_Gil"],
        y       = df_k[[metric]][df_k$site == "Rancho_San_Ignacio"],
        exact   = FALSE,
        correct = FALSE
      ),
      error = function(e) list(p.value = NA, statistic = NA)
    )
    stats_results[[paste(k, metric, sep = "_")]] <- data.frame(
      Kingdom = k,
      Metric  = metric,
      W       = ifelse(is.null(test_result$statistic), NA, test_result$statistic),
      p_value = test_result$p.value,
      Significativo = ifelse(!is.na(test_result$p.value) & test_result$p.value < 0.05, "*", "ns")
    )
  }
}

stats_df <- do.call(rbind, stats_results)
rownames(stats_df) <- NULL
cat("  Resultados Wilcoxon (Rancho Gil vs. Rancho San Ignacio):\n")
print(stats_df)

# ---- 5. Visualización: Boxplots de diversidad alfa ----
cat("\n===== GENERANDO FIGURAS =====\n")

# Función para crear boxplot de una métrica
plot_alpha_metric <- function(data, metric, title = NULL, y_label = NULL) {
  comparisons <- list(c("Rancho_Gil", "Rancho_San_Ignacio"))

  p <- ggplot(data, aes(x = site, y = .data[[metric]], fill = site, color = site)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
    scale_fill_manual(values  = site_colors, labels = site_labels) +
    scale_color_manual(values = site_colors, labels = site_labels) +
    scale_x_discrete(labels = site_labels) +
    labs(
      title = title %||% metric,
      x     = NULL,
      y     = y_label %||% metric,
      fill  = "Sitio",
      color = "Sitio"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position   = "none",
      axis.text.x       = element_text(angle = 25, hjust = 1, size = 10),
      plot.title        = element_text(face = "bold", size = 11),
      panel.grid.minor  = element_blank()
    ) +
    stat_compare_means(
      comparisons = comparisons,
      method      = "wilcox.test",
      label       = "p.format",
      size        = 3.5,
      tip.length  = 0.01
    )

  return(p)
}

# Crear plots individuales para bacterias
plots_bact <- lapply(c("Richness", "Shannon", "Simpson", "Evenness"), function(m) {
  plot_alpha_metric(
    data    = alpha_bact,
    metric  = m,
    title   = m,
    y_label = switch(m,
      Richness = "Riqueza observada (# taxa)",
      Shannon  = "Índice de Shannon (H')",
      Simpson  = "Índice de Simpson (1-D)",
      Evenness = "Equitatividad de Pielou (J)"
    )
  )
})

# Panel combinado: Bacterias
fig_alpha_bact <- wrap_plots(plots_bact, ncol = 4) +
  plot_annotation(
    title    = "Diversidad Alfa - Bacterias/Arqueas (16S rRNA, región V3-V4)",
    subtitle = "Suelo de vid (Vitis vinifera) | Baja California, México",
    caption  = "Prueba de Wilcoxon Mann-Whitney; n=3 por sitio",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      plot.caption  = element_text(size = 8,  color = "grey50")
    )
  )

ggsave(file.path(OUT_DIR, "alpha_diversity_bacteria.pdf"),
       fig_alpha_bact, width = 14, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "alpha_diversity_bacteria.png"),
       fig_alpha_bact, width = 14, height = 5, dpi = 300)
cat("  ✓ Figura guardada: alpha_diversity_bacteria.pdf/png\n")

# Panel combinado: Hongos
plots_fungi <- lapply(c("Richness", "Shannon", "Simpson", "Evenness"), function(m) {
  plot_alpha_metric(
    data    = alpha_fungi,
    metric  = m,
    title   = m,
    y_label = switch(m,
      Richness = "Riqueza observada (# taxa)",
      Shannon  = "Índice de Shannon (H')",
      Simpson  = "Índice de Simpson (1-D)",
      Evenness = "Equitatividad de Pielou (J)"
    )
  )
})

fig_alpha_fungi <- wrap_plots(plots_fungi, ncol = 4) +
  plot_annotation(
    title    = "Diversidad Alfa - Hongos (ITS1/ITS2)",
    subtitle = "Suelo de vid (Vitis vinifera) | Baja California, México",
    caption  = "Prueba de Wilcoxon Mann-Whitney; n=3 por sitio",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      plot.caption  = element_text(size = 8,  color = "grey50")
    )
  )

ggsave(file.path(OUT_DIR, "alpha_diversity_fungi.pdf"),
       fig_alpha_fungi, width = 14, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "alpha_diversity_fungi.png"),
       fig_alpha_fungi, width = 14, height = 5, dpi = 300)
cat("  ✓ Figura guardada: alpha_diversity_fungi.pdf/png\n")

# ---- 6. Curvas de rarefacción ----
cat("\n===== CURVAS DE RAREFACCIÓN =====\n")

plot_rarefaction <- function(pseudo_mat, meta_df, title_str, colors_vec) {
  # pseudo_mat: pseudo-conteos (taxa en filas, muestras en columnas)
  pseudo_t <- t(pseudo_mat)

  # Calcular curva de rarefacción con vegan
  # rarecurve() traza automáticamente, pero aquí extraemos los datos
  step <- max(1, floor(max(rowSums(pseudo_t)) / 50))  # 50 puntos en la curva
  rare_data <- rarecurve(pseudo_t, step = step, tidy = TRUE)
  
  rare_data <- rare_data |> dplyr::rename("Size"=Sample,"Sample" = Site)

  # rarecurve tidy=TRUE returns Sample as integer index; remap to sample names
  # sample_names <- rownames(pseudo_t)
  # rare_data$Sample <- sample_names[as.integer(rare_data$Sample)]

  # Añadir metadatos
  meta_df <- meta_df %>% rownames_to_column("Sample") %>% select(Sample, site)
  
  any(rare_data$Sample %in% meta_df$Sample)
  
  rare_data <- rare_data %>%
    left_join(
      meta_df,
      by = "Sample"
    )

  p <- ggplot(rare_data, aes(x = Size, y = Species,
                              group = Sample, color = site)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(data = rare_data %>% group_by(Sample) %>% slice_max(Sample, n = 1),
               size = 3) +
    scale_color_manual(values = colors_vec, labels = site_labels, name = "Sitio") +
    labs(
      title = title_str,
      x     = "Número de reads (pseudo-conteos)",
      y     = "Número de taxa observados"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11)
    )

  return(p)
}

# pseudo_mat, meta_df, title_str, colors_vec

fig_rare_bact <- plot_rarefaction(
  pseudo_mat = bact_pseudo,
  meta_df    = meta_clean,
  title_str  = "Curva de rarefacción - Bacterias/Arqueas (16S)",
  colors_vec = site_colors
)

fig_rare_fungi <- plot_rarefaction(
  pseudo_mat = fungi_pseudo,
  meta_df    = meta_clean,
  title_str  = "Curva de rarefacción - Hongos (ITS)",
  colors_vec = site_colors
)

fig_rarefaction <- fig_rare_bact / fig_rare_fungi +
  plot_annotation(
    caption = "Suelo de vid (Vitis vinifera) | Baja California, México",
    theme   = theme(plot.caption = element_text(color = "grey50"))
  )

ggsave(file.path(OUT_DIR, "rarefaction_curves.pdf"),
       fig_rarefaction, width = 10, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "rarefaction_curves.png"),
       fig_rarefaction, width = 10, height = 8, dpi = 300)
cat("  ✓ Figura guardada: rarefaction_curves.pdf/png\n")

# ---- 7. Exportar tabla de resultados ----
write.csv(alpha_all,  file = file.path(OUT_DIR, "alpha_diversity_all.csv"),  row.names = FALSE)
write.csv(stats_df,   file = file.path(OUT_DIR, "alpha_statistics.csv"),     row.names = FALSE)
cat("  ✓ Tablas exportadas: alpha_diversity_all.csv, alpha_statistics.csv\n")

# ---- 8. Resumen ----
cat("\n===== RESUMEN DIVERSIDAD ALFA =====\n")
cat("\nMedianas por sitio (Bacterias):\n")
alpha_bact %>%
  group_by(site) %>%
  summarise(across(c(Richness, Shannon, Simpson, Evenness, Chao1), median, na.rm = TRUE)) %>%
  print()

cat("\nMedianas por sitio (Hongos):\n")
alpha_fungi %>%
  group_by(site) %>%
  summarise(across(c(Richness, Shannon, Simpson, Evenness, Chao1), median, na.rm = TRUE)) %>%
  print()

cat("\n✓ Análisis de diversidad alfa completado. Procede con 03_beta_diversity.R\n")
