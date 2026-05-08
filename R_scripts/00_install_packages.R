# =============================================================================
# 00_install_packages.R
# Proyecto: Microbioma de Vid (Vitis vinifera) - Análisis de Amplicones 16S/ITS
# Autor: Ricardo Gómez-Reyes | rgomez41@uabc.edu.mx
# Fecha: 2026-04-26
#
# Descripción: Instalación de todos los paquetes R necesarios para el pipeline
#              de análisis de microbioma. Ejecutar UNA SOLA VEZ antes de correr
#              el resto de los scripts.
#
# Referencia: Wen et al. (2023) "The best practice for microbiome analysis
#             using R". Protein & Cell, 14, 713–725.
# =============================================================================

# ---- 1. Configurar repositorios ----
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- 2. Instalar BiocManager (gestor de paquetes Bioconductor) ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.22", ask = FALSE)  # Ajustar según versión de R

# ---- 3. Paquetes CRAN - Manipulación y lectura de datos ----
cran_packages <- c(
  # Lectura de datos
  "readxl",       # Leer archivos .xlsx
  "openxlsx",     # Leer y escribir .xlsx
  "readr",        # Lectura rápida de tablas

  # Manipulación de datos (tidyverse)
  "tidyverse",    # Colección: dplyr, ggplot2, tidyr, purrr, etc.
  "dplyr",        # Manipulación de data frames
  "tidyr",        # Transformación largo/ancho
  "reshape2",     # melt() y dcast() para reformatear datos
  "plyr",         # Split-Apply-Combine
  "stringr",      # Manipulación de cadenas de texto
  "forcats",      # Manejo de factores

  # Estadística multivariada y ecología
  "vegan",        # Análisis de comunidades: diversidad, ordenación, PERMANOVA
  "ape",          # Filogenética, PCoA, Mantel test
  "picante",      # Diversidad filogenética
  "ade4",         # Análisis multivariado, RDA, CCA
  "adegraphics",  # Visualización para ade4
  "adespatial",   # Análisis espacial y beta-diversidad
  "GUniFrac",     # UniFrac y rarefacción

  # Visualización
  "ggplot2",      # Sistema de gráficos principal
  "ggpubr",       # Gráficos de publicación (boxplots con stats)
  "pheatmap",     # Heatmaps
  "RColorBrewer", # Paletas de colores
  "viridis",      # Paletas de colores perceptualmente uniformes
  "cowplot",      # Combinación de múltiples gráficos
  "patchwork",    # Combinación modular de ggplots
  "ggrepel",      # Etiquetas que no se superponen en ggplot2
  "ggsci",        # Paletas de colores de revistas científicas
  "corrplot",     # Matrices de correlación
  "factoextra",   # Visualización de PCA/clustering
  "ggVennDiagram",# Diagramas de Venn con ggplot2
  "VennDiagram",  # Diagramas de Venn clásicos
  "scales",       # Escalas para ggplot2

  # Redes de co-ocurrencia
  "igraph",       # Análisis y visualización de redes
  "ggraph",       # Visualización de redes con ggplot2
  "tidygraph",    # Manipulación de grafos en estilo tidyverse
  "psych",        # Correlaciones y psicometría
  "Hmisc",        # Correlaciones con p-valores

  # Machine Learning y biomarkers
  "randomForest", # Random Forest
  "caret",        # Framework de machine learning
  "pROC",         # Curvas ROC
  "e1071",        # SVM y Naive Bayes

  # Utilidades
  "here",         # Rutas relativas portables
  "janitor",      # Limpieza de nombres de columnas
  "knitr",        # Reportes dinámicos
  "rmarkdown"     # Documentos R Markdown
)

# Instalar paquetes CRAN faltantes
installed_cran <- installed.packages()[,"Package"]
to_install_cran <- cran_packages[!cran_packages %in% installed_cran]

if (length(to_install_cran) > 0) {
  cat("Instalando", length(to_install_cran), "paquetes de CRAN...\n")
  install.packages(to_install_cran, dependencies = TRUE)
} else {
  cat("Todos los paquetes CRAN ya están instalados.\n")
}

# ---- 4. Paquetes Bioconductor ----
bioc_packages <- c(
  "phyloseq",         # Objeto integrador de datos de microbioma (FUNDAMENTAL)
  "microbiome",       # Extensión de phyloseq con funciones adicionales
  "DESeq2",           # Abundancia diferencial (modelos binomiales negativos)
  "edgeR",            # Abundancia diferencial (TMM normalization)
  "limma",            # Abundancia diferencial (modelos lineales)
  "ALDEx2",           # Abundancia diferencial (mejores resultados en benchmarks)
  "metagenomeSeq",    # Normalización CSS y abundancia diferencial
  "Biostrings",       # Manejo de secuencias biológicas
  "biomformat",       # Lectura de archivos .biom (QIIME)
  "ggtree",           # Visualización de árboles filogenéticos
  "treeio",           # Lectura de árboles filogenéticos
  "ANCOMBC",          # Análisis de composición con corrección de bias
  "MicrobiomeAnalystR" # Suite de análisis de microbioma (instalar si está disponible)
)

installed_bioc <- installed.packages()[,"Package"]
to_install_bioc <- bioc_packages[!bioc_packages %in% installed_bioc]

if (length(to_install_bioc) > 0) {
  cat("Instalando", length(to_install_bioc), "paquetes de Bioconductor...\n")
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
} else {
  cat("Todos los paquetes Bioconductor ya están instalados.\n")
}

# ---- 5. Paquetes de GitHub (opcionales, requieren devtools/remotes) ----
# Descomenta los que necesites

# install.packages("remotes")
# remotes::install_github("taowenmicro/EasyMicrobiomeR")   # Pipeline integrado del paper
# remotes::install_github("zdk123/SpiecEasi")               # Redes de co-ocurrencia (SparCC)
# remotes::install_github("microsud/microbiomeutilities")  # Utilidades adicionales
# remotes::install_github("vmikk/metagMisc")               # Utilidades metagenómicas

# ---- 6. Verificación de instalación ----

cat("\n========== VERIFICACIÓN DE PAQUETES ==========\n")
all_packages <- c(cran_packages, bioc_packages)
installed_all <- installed.packages()[,"Package"]

missing <- all_packages[!all_packages %in% installed_all]
installed_ok <- all_packages[all_packages %in% installed_all]

cat("✓ Instalados:", length(installed_ok), "/", length(all_packages), "\n")
if (length(missing) > 0) {
  cat("✗ Faltantes:", paste(missing, collapse = ", "), "\n")
  cat("  → Intenta instalarlos manualmente con:\n")
  cat("    install.packages(c('", paste(missing, collapse = "', '"), "'))\n")
} else {
  cat("✓ ¡Todos los paquetes instalados correctamente!\n")
}

# ---- 7. Cargar paquetes clave y verificar versiones ----
cat("\n========== VERSIONES CLAVE ==========\n")
key_packages <- c("vegan", "phyloseq", "ggplot2", "dplyr", "DESeq2", "ALDEx2")
for (pkg in key_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-15s v%s\n", pkg, packageVersion(pkg)))
  }
}

cat("\n✓ Setup completado. Procede con 01_data_import.R\n")
