# =============================================================================
# 01_data_import.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Importación, limpieza, mapeo de muestras y creación de objetos
#              phyloseq para bacterias/arqueas (16S) y hongos (ITS).
#              Este script es la BASE de todo el pipeline - los demás scripts
#              lo llaman con source("01_data_import.R").
#
# Datos de entrada:
#   - bact_abund.xlsx  : Abundancias relativas bacterias/arqueas (827 taxa × 6 muestras)
#   - fungi_abund.xlsx : Abundancias relativas hongos (582 taxa × 6 muestras)
#   - datos_suelo.xlsx : Metadatos fisicoquímicos del suelo (6 muestras × 21 variables)
#
# NOTA IMPORTANTE SOBRE MUESTRAS:
#   Las matrices de abundancia usan IDs: DRI000–DRI005
#   Los metadatos usan nombres:          RANCHO_GIL_1–3, RANCHO_SAN_IGNACIO_1–3
#   El mapeo asumido es:
#     DRI000 = RANCHO_GIL_1       (Rancho Gil, sitio 1)
#     DRI001 = RANCHO_GIL_2       (Rancho Gil, sitio 2)
#     DRI002 = RANCHO_GIL_3       (Rancho Gil, sitio 3)
#     DRI003 = RANCHO_SAN_IGNACIO_1 (Rancho San Ignacio, sitio 1)
#     DRI004 = RANCHO_SAN_IGNACIO_2 (Rancho San Ignacio, sitio 2)
#     DRI005 = RANCHO_SAN_IGNACIO_3 (Rancho San Ignacio, sitio 3)
#   → VERIFICAR este mapeo con el laboratorio antes del análisis final.
# =============================================================================


rm(list = ls())

if(!is.null(dev.list())) dev.off()

# ---- 0. Paquetes requeridos ----
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(phyloseq)
})

# ---- 1. Definir rutas ----
# Adaptar la ruta base al directorio de trabajo del proyecto
DATA_DIR <- here::here("~/Documents/Metagenomics/Microbial-soil/Artículo/")  # Ajustar si los scripts están en R_scripts/

# Alternativa con rutas absolutas (descomentar si es necesario):
# DATA_DIR <- "/Users/rjegr/Documents/Metagenomics/Microbial-soil/Artículo"

# Archivos de entrada
FILE_BACT  <- file.path(DATA_DIR, "bact_abund.xlsx")
FILE_FUNGI <- file.path(DATA_DIR, "fungi_abund.xlsx")
FILE_META  <- file.path(DATA_DIR, "datos_suelo.xlsx")

# ---- 2. Mapeo de IDs de muestra ----
# IMPORTANTE: Verificar con el proveedor de secuenciación (BeCrop / laboratorio)
sample_map <- data.frame(
  sample_id  = c("DRI000", "DRI001", "DRI002", "DRI003", "DRI004", "DRI005"),
  site_name  = c("RANCHO_GIL_1", "RANCHO_GIL_2", "RANCHO_GIL_3",
                 "RANCHO_SAN_IGNACIO_1", "RANCHO_SAN_IGNACIO_2", "RANCHO_SAN_IGNACIO_3"),
  site       = c("Rancho_Gil", "Rancho_Gil", "Rancho_Gil",
                 "Rancho_San_Ignacio", "Rancho_San_Ignacio", "Rancho_San_Ignacio"),
  replicate  = c(1, 2, 3, 1, 2, 3),
  stringsAsFactors = FALSE
)

# ---- 3. Función auxiliar: cargar matriz de abundancia ----
load_abund_matrix <- function(filepath, sheet = "abund") {
  # Leer hoja de Excel
  raw <- read_excel(filepath, sheet = sheet, col_names = TRUE)
  
  # La primera columna es el nombre del taxón
  taxa_col <- names(raw)[1]
  cat(sprintf("  → Columna de taxa: '%s'\n", taxa_col))

  # Eliminar filas con taxa = NA (filas vacías al final)
  raw <- raw %>% filter(!is.na(.data[[taxa_col]]))

  # Convertir a data.frame con taxa como rownames
  mat <- raw %>%
    column_to_rownames(var = taxa_col) %>%
    as.matrix()

  # Verificar que todas las columnas son numéricas
  if (!is.numeric(mat)) {
    rn <- rownames(mat)
    mat <- apply(mat, 2, as.numeric)
    rownames(mat) <- rn
  }

  # Reemplazar NA con 0
  mat[is.na(mat)] <- 0

  cat(sprintf("  → Dimensiones: %d taxa × %d muestras\n", nrow(mat), ncol(mat)))
  cat(sprintf("  → Rango de valores: [%.4f, %.4f]\n", min(mat), max(mat)))
  cat(sprintf("  → Suma total por muestra:\n"))
  print(round(colSums(mat), 2))

  return(mat)
}

# ---- 4. Función: extraer taxonomía desde nombres de taxa ----
# Los nombres son a nivel de especie (2 palabras = Género Especie)
# o con "sp." para morfoespecies no identificadas
build_tax_table <- function(taxa_names) {
  # Extraer género (primera palabra)
  genus <- str_extract(taxa_names, "^[A-Za-z]+")

  # Detectar si es morfoespecies (contiene "sp." o "cf.")
  is_morph <- str_detect(taxa_names, "\\bsp\\.?$|\\bcf\\.\\b|\\baff\\.\\b")

  # Crear tabla taxonómica
  tax_df <- data.frame(
    Genus   = genus,
    Species = taxa_names,
    is_morphospecies = is_morph,
    row.names = taxa_names,
    stringsAsFactors = FALSE
  )

  # Nota: la tabla de taxonomía completa (Phylum, Class, Order, Family)
  # debe obtenerse de la base de datos SILVA (16S) o UNITE (ITS)
  # mediante QIIME2 o DADA2. Esta es una versión simplificada.

  return(tax_df)
}

# ---- 5. Función: convertir a pseudo-conteos (para métodos que requieren counts) ----
# Como los datos son abundancias relativas (%), multiplicamos para obtener
# pseudo-conteos enteros. Factor = 10000 simula ~10,000 reads por muestra.
relabund_to_pseudocounts <- function(mat, factor = 10000) {
  pseudo <- round(mat * factor / 100)
  # Asegurar mínimo de 0
  pseudo[pseudo < 0] <- 0
  return(pseudo)
}

# ---- 6. Cargar datos de bacterias/arqueas (16S) ----
cat("\n===== CARGANDO DATOS BACTERIAS/ARQUEAS (16S) =====\n")
bact_mat <- load_abund_matrix(FILE_BACT, sheet = "abund")

# Renombrar columnas con IDs cortos (DRI000-DRI005)
# Verificar que los nombres coinciden con sample_map
colnames(bact_mat)  # Debe ser: DRI000 DRI001 DRI002 DRI003 DRI004 DRI005

# ---- 7. Cargar datos de hongos (ITS) ----
cat("\n===== CARGANDO DATOS HONGOS (ITS) =====\n")
fungi_mat <- load_abund_matrix(FILE_FUNGI, sheet = "abund")

# ---- 8. Cargar metadatos de suelo ----
cat("\n===== CARGANDO METADATOS DE SUELO =====\n")
meta_raw <- read_excel(FILE_META, sheet = "fis_qui")
cat(sprintf("  → Dimensiones: %d muestras × %d variables\n", nrow(meta_raw), ncol(meta_raw)))
cat("  → Variables:\n")
cat(paste("   ", names(meta_raw), collapse = "\n"), "\n")

# Limpiar nombres de columnas (reemplazar espacios y caracteres especiales)
meta_clean <- meta_raw %>%
  rename(site_name = Muestra) %>%
  # Unir con el mapeo de IDs
  left_join(sample_map, by = "site_name") %>%
  # sample_id como identificador de fila
  column_to_rownames(var = "sample_id")

# Renombrar columnas para facilitar el análisis
# (nombres sin unidades para uso en funciones de R)
names(meta_clean) <- names(meta_clean) %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("[^A-Za-z0-9_]", "") %>%
  str_to_lower()

cat("  → Columnas limpias:\n")
print(names(meta_clean))
print(head(meta_clean[, 1:8]))

# ---- 9. Crear tabla taxonómica ----
cat("\n===== CONSTRUYENDO TABLAS TAXONÓMICAS =====\n")
tax_bact  <- build_tax_table(rownames(bact_mat))
tax_fungi <- build_tax_table(rownames(fungi_mat))
cat(sprintf("  → Taxa bacterias: %d | Taxa hongos: %d\n", nrow(tax_bact), nrow(tax_fungi)))

# ---- 10. Crear objetos phyloseq ----
cat("\n===== CREANDO OBJETOS PHYLOSEQ =====\n")

# Asegurarse de que el orden de muestras coincida entre matriz y metadata
# (los colnames de la matriz deben ser los rownames de metadata)
samples_in_meta <- rownames(meta_clean)
samples_in_bact <- colnames(bact_mat)

if (!all(samples_in_bact %in% samples_in_meta)) {
  warning("¡ATENCIÓN! Algunos samples de la matriz de bacterias no están en los metadatos.")
  cat("  En bact pero no en meta:", setdiff(samples_in_bact, samples_in_meta), "\n")
}

# Reordenar metadata para que coincida con el orden de columnas de la matriz
meta_ordered_bact  <- meta_clean[colnames(bact_mat), , drop = FALSE]
meta_ordered_fungi <- meta_clean[colnames(fungi_mat), , drop = FALSE]

# Crear pseudo-conteos (para métodos que requieren counts enteros)
bact_pseudo  <- relabund_to_pseudocounts(bact_mat,  factor = 10000)
fungi_pseudo <- relabund_to_pseudocounts(fungi_mat, factor = 10000)

# --- Objeto phyloseq con ABUNDANCIAS RELATIVAS (%) ---
ps_bact <- phyloseq(
  otu_table(bact_mat,         taxa_are_rows = TRUE),
  sample_data(meta_ordered_bact),
  tax_table(as.matrix(tax_bact))
)

ps_fungi <- phyloseq(
  otu_table(fungi_mat,        taxa_are_rows = TRUE),
  sample_data(meta_ordered_fungi),
  tax_table(as.matrix(tax_fungi))
)

# --- Objeto phyloseq con PSEUDO-CONTEOS (para DESeq2, etc.) ---
ps_bact_counts <- phyloseq(
  otu_table(bact_pseudo,      taxa_are_rows = TRUE),
  sample_data(meta_ordered_bact),
  tax_table(as.matrix(tax_bact))
)

ps_fungi_counts <- phyloseq(
  otu_table(fungi_pseudo,     taxa_are_rows = TRUE),
  sample_data(meta_ordered_fungi),
  tax_table(as.matrix(tax_fungi))
)

cat("  ✓ phyloseq bacterias (rel. abund.):\n"); print(ps_bact)
cat("  ✓ phyloseq hongos (rel. abund.):\n");    print(ps_fungi)

# ---- 11. Filtrado básico de taxa ----
# Eliminar taxa con abundancia 0 en todas las muestras (taxa fantasma)
cat("\n===== FILTRADO BÁSICO =====\n")

ps_bact  <- prune_taxa(taxa_sums(ps_bact)  > 0, ps_bact)
ps_fungi <- prune_taxa(taxa_sums(ps_fungi) > 0, ps_fungi)
ps_bact_counts  <- prune_taxa(taxa_sums(ps_bact_counts)  > 0, ps_bact_counts)
ps_fungi_counts <- prune_taxa(taxa_sums(ps_fungi_counts) > 0, ps_fungi_counts)

cat(sprintf("  ✓ Bacterias tras filtrado: %d taxa\n", ntaxa(ps_bact)))
cat(sprintf("  ✓ Hongos tras filtrado:    %d taxa\n", ntaxa(ps_fungi)))

# Prevalencia: taxa presentes en al menos 2 muestras (33%)
prev_threshold <- 2
ps_bact_prev  <- filter_taxa(ps_bact,  function(x) sum(x > 0) >= prev_threshold, TRUE)
ps_fungi_prev <- filter_taxa(ps_fungi, function(x) sum(x > 0) >= prev_threshold, TRUE)

cat(sprintf("  ✓ Bacterias (prevalencia ≥ %d muestras): %d taxa\n", prev_threshold, ntaxa(ps_bact_prev)))
cat(sprintf("  ✓ Hongos    (prevalencia ≥ %d muestras): %d taxa\n", prev_threshold, ntaxa(ps_fungi_prev)))

# ---- 12. Resumen estadístico básico ----
cat("\n===== RESUMEN DE LOS DATOS =====\n")
cat("\n-- Metadatos (primeras filas) --\n")
print(meta_clean %>% select(site, replicate, everything()) %>% head())

cat("\n-- Distribución de muestras por sitio --\n")
print(table(meta_clean$site))

cat("\n-- Riqueza observada (taxa con abund > 0) --\n")
bact_richness  <- colSums(bact_mat > 0)
fungi_richness <- colSums(fungi_mat > 0)
cat("  Bacterias:\n"); print(bact_richness)
cat("  Hongos:\n");    print(fungi_richness)

# ---- 13. Exportar objetos para uso en otros scripts ----
# Los objetos creados en este script están disponibles en el entorno global
# cuando se usa: source("01_data_import.R")

cat("\n===== OBJETOS DISPONIBLES PARA ANÁLISIS =====\n")
cat("
  Objetos cargados en el entorno R:
  ─────────────────────────────────────────────────────────────────
  sample_map      : Mapeo DRI_ID ↔ Rancho (data.frame, 6 × 4)
  bact_mat        : Matriz de abundancias bacterianas en % (827 × 6)
  fungi_mat       : Matriz de abundancias fúngicas en % (582 × 6)
  bact_pseudo     : Pseudo-conteos bacterias (× 10,000) (827 × 6)
  fungi_pseudo    : Pseudo-conteos hongos (× 10,000) (582 × 6)
  meta_clean      : Metadatos fisicoquímicos limpios (6 × 21+)
  tax_bact        : Tabla taxonómica bacterias (Genus + Species)
  tax_fungi       : Tabla taxonómica hongos (Genus + Species)

  Objetos phyloseq:
  ps_bact         : Bacterias, abundancias relativas %
  ps_fungi        : Hongos, abundancias relativas %
  ps_bact_counts  : Bacterias, pseudo-conteos (para DESeq2)
  ps_fungi_counts : Hongos, pseudo-conteos (para DESeq2)
  ps_bact_prev    : Bacterias filtradas por prevalencia (≥ 2 muestras)
  ps_fungi_prev   : Hongos filtrados por prevalencia (≥ 2 muestras)
  ─────────────────────────────────────────────────────────────────
  NOTA: n=6 muestras (3 Rancho Gil + 3 Rancho San Ignacio).
        El poder estadístico es limitado; interpretar con precaución.
")

# ---- 14. Guardar workspace para uso offline ----
# (Descomentar para guardar todos los objetos en un archivo .RData)
# save(list = c("ps_bact", "ps_fungi", "ps_bact_counts", "ps_fungi_counts",
#               "ps_bact_prev", "ps_fungi_prev",
#               "bact_mat", "fungi_mat", "meta_clean", "sample_map",
#               "tax_bact", "tax_fungi"),
#      file = "microbioma_vid_data.RData")
# cat("  → Datos guardados en 'microbioma_vid_data.RData'\n")

cat("\n✓ data_import completado. Procede con 02_alpha_diversity.R\n")
