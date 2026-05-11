# =============================================================================
# 07_network_analysis.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Análisis de redes de co-ocurrencia microbiana (co-occurrence
#              networks) para identificar patrones de asociación entre taxa.
#
# Análisis incluidos:
#   1. Cálculo de correlaciones (Spearman, con corrección múltiple)
#   2. Construcción de redes de co-ocurrencia (igraph)
#   3. Visualización de redes (ggraph)
#   4. Métricas topológicas de la red (grado, betweenness, clustering)
#   5. Detección de módulos/comunidades en la red (Louvain)
#   6. Identificación de nodos hub (keystone taxa)
#   7. Red inter-reino (bacterias + hongos)
#
# NOTA IMPORTANTE: Las redes de co-ocurrencia son más robustas con n≥10 muestras.
#   Con n=6 muestras, estas redes son ALTAMENTE EXPLORATORIAS y deben
#   interpretarse con precaución. Todos los resultados son tentativas.
#   Se recomienda esperar más réplicas antes de publicar.
#
# Dependencias: 01_data_import.R
# =============================================================================

rm(list = ls())

if(!is.null(dev.list())) dev.off()

# ---- 0. Cargar datos ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("01_data_import.R")

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(psych)      # corr.test() con corrección de múltiples pruebas
  library(Hmisc)      # rcorr() rápido
  library(stringr)
  library(scales)
})

# ---- 1. Configuración ----

OUT_DIR1 <- here::here("~/Documents/Metagenomics/Microbial-soil/Artículo/")  # Ajustar si los scripts están en R_scripts/

OUT_DIR2 <- "figures/network_analysis"

OUT_DIR <- file.path(OUT_DIR1, OUT_DIR2)


dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Umbrales de red
COR_THRESH   <- 0.7    # Correlación mínima para un enlace (|r| ≥ 0.7)
P_THRESH     <- 0.05   # p-valor máximo (sin corrección, dado n pequeño)
MIN_PREV_NET <- 3      # Prevalencia mínima para incluir en la red (≥ 3/6 muestras)

# ---- 2. Preparar matrices para correlación ----
# Filtrar taxa por prevalencia mínima (reduce ruido)
filter_for_network <- function(mat, min_prev = 3, top_n = NULL) {
  # Filtrar por prevalencia
  keep <- rowSums(mat > 0) >= min_prev
  mat_f <- mat[keep, , drop = FALSE]
  cat(sprintf("  Prevalencia ≥ %d: %d / %d taxa\n", min_prev, sum(keep), nrow(mat)))

  # Opcionalmente, conservar solo los top-N por abundancia media
  if (!is.null(top_n) && nrow(mat_f) > top_n) {
    top_idx <- order(rowMeans(mat_f), decreasing = TRUE)[1:top_n]
    mat_f   <- mat_f[top_idx, , drop = FALSE]
    cat(sprintf("  Reducido a Top %d por abundancia media\n", top_n))
  }

  return(mat_f)
}

bact_net  <- filter_for_network(bact_mat,  min_prev = MIN_PREV_NET, top_n = 100)
fungi_net <- filter_for_network(fungi_mat, min_prev = MIN_PREV_NET, top_n = 100)

# ---- 3. Función: calcular correlaciones Spearman ----
compute_correlations <- function(mat, method = "spearman") {
  # mat: taxa en filas, muestras en columnas
  # Transponer: muestras en filas, taxa en columnas
  mat_t <- t(mat)

  # Usar psych::corr.test() para obtener r y p-valores
  # NOTA: con n=6, los p-valores tienen muy poco poder
  cor_result <- corr.test(mat_t, method = method, adjust = "BH", ci = FALSE)

  return(list(
    r    = cor_result$r,     # Matriz de correlaciones
    p    = cor_result$p,     # p-valores ajustados (BH)
    p_raw = cor_result$p.adj # p-valores sin ajustar (para exploración)
  ))
}

cat("\n===== CALCULANDO CORRELACIONES =====\n")
cat("  Bacterias...\n")
cor_bact  <- compute_correlations(bact_net)
cat("  Hongos...\n")
cor_fungi <- compute_correlations(fungi_net)

# ---- 4. Función: construir grafo a partir de correlaciones ----
build_network <- function(cor_mat, p_mat, taxa_names,
                           r_thresh = COR_THRESH, p_thresh = P_THRESH,
                           kingdom = "Bacteria") {
  n <- nrow(cor_mat)

  # Crear lista de aristas (edges) que pasan el umbral
  edges <- data.frame()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      r_val <- cor_mat[i, j]
      p_val <- p_mat[i, j]

      if (!is.na(r_val) && !is.na(p_val) &&
          abs(r_val) >= r_thresh && p_val <= p_thresh) {
        edges <- rbind(edges, data.frame(
          from      = taxa_names[i],
          to        = taxa_names[j],
          weight    = abs(r_val),
          r         = r_val,
          p_adj     = p_val,
          type      = ifelse(r_val > 0, "positive", "negative"),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(edges) == 0) {
    cat(sprintf("  ⚠ Sin aristas con |r| ≥ %.1f y p ≤ %.2f. Reduciendo umbrales...\n",
                r_thresh, p_thresh))
    # Umbral más permisivo
    r_thresh <- 0.5
    return(build_network(cor_mat, p_mat, taxa_names, r_thresh, p_thresh * 2, kingdom))
  }

  cat(sprintf("  ✓ Red %s: %d nodos, %d aristas (|r| ≥ %.1f, p ≤ %.2f)\n",
              kingdom, length(unique(c(edges$from, edges$to))),
              nrow(edges), r_thresh, p_thresh))

  # Crear grafo
  g <- graph_from_data_frame(edges, directed = FALSE)

  # Añadir atributos de nodo
  V(g)$kingdom <- kingdom
  V(g)$genus   <- str_extract(V(g)$name, "^[A-Za-z]+")

  return(g)
}

cat("\n===== CONSTRUYENDO REDES =====\n")
g_bact  <- build_network(cor_bact$r,  cor_bact$p,  rownames(cor_bact$r),  kingdom = "Bacteria/Archaea")
g_fungi <- build_network(cor_fungi$r, cor_fungi$p, rownames(cor_fungi$r), kingdom = "Fungi")

# ---- 5. Función: calcular métricas topológicas ----
compute_network_metrics <- function(g, label = "") {
  # Métricas de nodo
  metrics <- data.frame(
    node          = V(g)$name,
    kingdom       = V(g)$kingdom,
    genus         = V(g)$genus,
    degree        = degree(g),
    betweenness   = round(betweenness(g, normalized = TRUE), 4),
    closeness     = round(closeness(g, normalized = TRUE), 4),
    strength      = round(strength(g), 4),  # suma de pesos de aristas
    stringsAsFactors = FALSE
  )

  # Módulos (Louvain community detection)
  set.seed(42)
  community   <- cluster_louvain(g)
  metrics$module <- membership(community)

  # Clasificar como hub si degree ≥ percentil 75
  deg_75 <- quantile(metrics$degree, 0.75)
  metrics$is_hub <- metrics$degree >= deg_75

  # Métricas globales
  global <- list(
    label               = label,
    n_nodes             = vcount(g),
    n_edges             = ecount(g),
    n_positive          = sum(E(g)$r > 0),
    n_negative          = sum(E(g)$r < 0),
    density             = round(edge_density(g), 4),
    avg_degree          = round(mean(degree(g)), 2),
    avg_clustering      = round(transitivity(g, type = "average"), 4),
    modularity          = round(modularity(community), 4),
    n_modules           = max(membership(community)),
    diameter            = diameter(g),
    avg_path_length     = round(mean_distance(g), 3)
  )

  cat(sprintf("\n  === Métricas globales - %s ===\n", label))
  cat(sprintf("  Nodos: %d | Aristas: %d (+ %d | - %d)\n",
              global$n_nodes, global$n_edges, global$n_positive, global$n_negative))
  cat(sprintf("  Densidad: %.4f | Grado medio: %.2f\n",
              global$density, global$avg_degree))
  cat(sprintf("  Coef. clustering: %.4f | Modularidad: %.4f (%d módulos)\n",
              global$avg_clustering, global$modularity, global$n_modules))
  cat(sprintf("  Diámetro: %.0f | Distancia media: %.3f\n",
              global$diameter, global$avg_path_length))

  cat(sprintf("\n  Top 10 hubs (por grado):\n"))
  print(head(metrics %>% arrange(desc(degree)), 10) %>%
          select(node, degree, betweenness, module))

  return(list(node_metrics = metrics, global_metrics = global))
}

cat("\n===== MÉTRICAS DE RED =====\n")
metrics_bact  <- compute_network_metrics(g_bact,  "Bacterias/Arqueas")
metrics_fungi <- compute_network_metrics(g_fungi, "Hongos")

# ---- 6. Función: visualizar red ----
plot_network <- function(g, node_metrics, title_str, filename,
                          max_label_degree = 5) {
  # Convertir a tidygraph para ggraph
  tg <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    left_join(node_metrics %>% rename(name = node), by = "name") %>%
    mutate(
      module_f = as.factor(module),
      label    = ifelse(degree >= sort(degree, decreasing = TRUE)[min(max_label_degree, n())],
                        str_trunc(name, 20, "right"), NA)
    ) %>%
    activate(edges) %>%
    mutate(edge_color = ifelse(r > 0, "positive", "negative"))

  # Paleta de módulos
  n_modules <- max(node_metrics$module, na.rm = TRUE)
  module_colors <- colorRampPalette(brewer.pal(min(n_modules, 9), "Set1"))(n_modules)

  p <- ggraph(tg, layout = "fr") +  # Fruchterman-Reingold layout
    # Aristas
    geom_edge_link(
      aes(color = edge_color, width = weight, alpha = weight),
      show.legend = TRUE
    ) +
    scale_edge_color_manual(
      values = c("positive" = "#e53935", "negative" = "#1e88e5"),
      labels = c("Co-ocurrencia (+)", "Exclusión mutua (-)"),
      name   = "Asociación"
    ) +
    scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    scale_edge_alpha(range = c(0.3, 0.8), guide = "none") +
    # Nodos
    geom_node_point(
      aes(size = degree, fill = module_f),
      shape = 21, color = "grey30", alpha = 0.85
    ) +
    geom_node_text(
      aes(label = label),
      size = 2, repel = TRUE,
      max.overlaps = 15, color = "grey10"
    ) +
    scale_size_continuous(range = c(2, 10), name = "Grado") +
    scale_fill_manual(values = module_colors, name = "Módulo") +
    labs(
      title    = title_str,
      subtitle = sprintf("n=%d nodos | %d aristas | %d módulos",
                         vcount(g), ecount(g),
                         max(node_metrics$module, na.rm = TRUE))
    ) +
    theme_graph(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "right"
    )

  ggsave(file.path(OUT_DIR, filename), p, width = 12, height = 9, dpi = 300, )
  cat(sprintf("  ✓ Red guardada: %s\n", filename))
  return(invisible(p))
}

cat("\n===== VISUALIZANDO REDES =====\n")
plot_network(g_bact,  metrics_bact$node_metrics,
             "Red de Co-ocurrencia — Bacterias/Arqueas (16S)",
             "network_bacteria.pdf")
plot_network(g_fungi, metrics_fungi$node_metrics,
             "Red de Co-ocurrencia — Hongos (ITS)",
             "network_fungi.pdf")

# ---- 7. Red inter-reino (bacterias + hongos) ----
cat("\n===== RED INTER-REINO (Bacterias × Hongos) =====\n")
cat("  Construyendo red de co-ocurrencia entre bacterias y hongos...\n")

# Calcular correlaciones cruzadas
cor_cross <- corr.test(t(bact_net), t(fungi_net), method = "spearman",
                        adjust = "BH", ci = FALSE)

# Construir aristas entre bacterias y hongos
edges_cross <- data.frame()
for (i in 1:nrow(cor_cross$r)) {
  for (j in 1:ncol(cor_cross$r)) {
    r_val <- cor_cross$r[i, j]
    p_val <- cor_cross$p[i, j]
    if (!is.na(r_val) && !is.na(p_val) &&
        abs(r_val) >= COR_THRESH && p_val <= P_THRESH) {
      edges_cross <- rbind(edges_cross, data.frame(
        from    = rownames(cor_cross$r)[i],
        to      = colnames(cor_cross$r)[j],
        weight  = abs(r_val),
        r       = r_val,
        p_adj   = p_val,
        type    = ifelse(r_val > 0, "positive", "negative"),
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat(sprintf("  Inter-reino: %d asociaciones inter-reino detectadas\n", nrow(edges_cross)))

if (nrow(edges_cross) > 0) {
  # Construir grafo inter-reino
  nodes_bact_cross  <- data.frame(name = rownames(bact_net), kingdom = "Bacteria/Archaea")
  nodes_fungi_cross <- data.frame(name = rownames(fungi_net), kingdom = "Fungi")
  nodes_all         <- rbind(nodes_bact_cross, nodes_fungi_cross)
  nodes_all         <- nodes_all[!duplicated(nodes_all$name), ]

  g_cross <- graph_from_data_frame(edges_cross, directed = FALSE,
                                    vertices = nodes_all[nodes_all$name %in%
                                                         c(edges_cross$from, edges_cross$to), ])

  tg_cross <- as_tbl_graph(g_cross) %>%
    activate(nodes) %>%
    mutate(degree = centrality_degree())

  p_cross <- ggraph(tg_cross, layout = "fr") +
    geom_edge_link(aes(color = type, width = weight, alpha = weight)) +
    scale_edge_color_manual(
      values = c("positive" = "#e53935", "negative" = "#1e88e5"),
      name   = "Asociación"
    ) +
    scale_edge_width(range = c(0.3, 1.5), guide = "none") +
    scale_edge_alpha(range = c(0.4, 0.9), guide = "none") +
    geom_node_point(aes(size = degree, fill = kingdom),
                    shape = 21, color = "grey30") +
    geom_node_text(aes(label = str_trunc(name, 18, "right")),
                   size = 1.8, repel = TRUE) +
    scale_fill_manual(
      values = c("Bacteria/Archaea" = "#2196F3", "Fungi" = "#FF9800"),
      name   = "Reino"
    ) +
    scale_size_continuous(range = c(2, 8), name = "Grado") +
    labs(
      title    = "Red Inter-reino — Co-ocurrencia Bacterias × Hongos",
      subtitle = sprintf("%d asociaciones entre bacterias y hongos detectadas",
                         ecount(g_cross))
    ) +
    theme_graph(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))

  ggsave(file.path(OUT_DIR, "network_inter_kingdom.pdf"),
         p_cross, width = 12, height = 9, dpi = 300, )
  cat("  ✓ Red inter-reino guardada: network_inter_kingdom.pdf\n")

  # Exportar aristas inter-reino
  write.csv(edges_cross,
            file.path(OUT_DIR, "inter_kingdom_associations.csv"), row.names = FALSE)
}

# ---- 8. Exportar resultados ----
cat("\n===== EXPORTANDO RESULTADOS =====\n")

write.csv(metrics_bact$node_metrics,
          file.path(OUT_DIR, "network_node_metrics_bacteria.csv"), row.names = FALSE)
write.csv(metrics_fungi$node_metrics,
          file.path(OUT_DIR, "network_node_metrics_fungi.csv"),    row.names = FALSE)

# Tabla de aristas
bact_edges <- as_data_frame(g_bact,  what = "edges")
fungi_edges <- as_data_frame(g_fungi, what = "edges")
write.csv(bact_edges,  file.path(OUT_DIR, "network_edges_bacteria.csv"), row.names = FALSE)
write.csv(fungi_edges, file.path(OUT_DIR, "network_edges_fungi.csv"),    row.names = FALSE)

# Resumen global
global_summary <- data.frame(
  Red = c("Bacterias/Arqueas", "Hongos"),
  Nodos    = c(metrics_bact$global$n_nodes,  metrics_fungi$global$n_nodes),
  Aristas  = c(metrics_bact$global$n_edges,  metrics_fungi$global$n_edges),
  Positivas = c(metrics_bact$global$n_positive, metrics_fungi$global$n_positive),
  Negativas = c(metrics_bact$global$n_negative, metrics_fungi$global$n_negative),
  Densidad  = c(metrics_bact$global$density,  metrics_fungi$global$density),
  Grado_medio = c(metrics_bact$global$avg_degree, metrics_fungi$global$avg_degree),
  Clustering  = c(metrics_bact$global$avg_clustering, metrics_fungi$global$avg_clustering),
  Modularidad = c(metrics_bact$global$modularity, metrics_fungi$global$modularity),
  N_Modulos   = c(metrics_bact$global$n_modules,  metrics_fungi$global$n_modules)
)
write.csv(global_summary, file.path(OUT_DIR, "network_global_metrics.csv"), row.names = FALSE)
cat("  ✓ Métricas globales exportadas: network_global_metrics.csv\n")
print(global_summary)

cat("\n
  ⚠ ADVERTENCIA FINAL - REDES DE CO-OCURRENCIA:
  ─────────────────────────────────────────────────────────────────
  Con n=6 muestras, las correlaciones de Spearman tienen muy bajo
  poder estadístico y alta tasa de falsos positivos. Las redes
  aquí generadas son EXPLORATORIAS.

  Para análisis robustos se recomienda:
  • n ≥ 20 muestras (idealmente n ≥ 50)
  • Usar SparCC (SpiecEasi) en lugar de Spearman para datos composicionales
  • Bootstrap o permutaciones para validar aristas
  ─────────────────────────────────────────────────────────────────
")

cat("\n✓ Análisis de redes de co-ocurrencia completado.\n")
cat("  Pipeline completo finalizado. Revisa la carpeta 'figures/' para todos los resultados.\n")
