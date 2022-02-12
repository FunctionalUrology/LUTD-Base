# ########################################################################
#                             Setting up                                 #
# ########################################################################

set.seed(45)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
options(connectionObserver = NULL)

#load and install github packages separately
#devtools::install_github("kassambara/easyGgplot2")
#devtools::install_github("ssayols/rrvgo")
#devtools::install_github("jokergoo/InteractiveComplexHeatmap")
#devtools::install_github("jokergoo/ComplexHeatmap")

packages <- c("pheatmap","devtools","biomaRt","ggforce","tidyverse","easyGgplot2",
              "clusterProfiler", "circlize","dplyr","DESeq2","plotly","ggridges",
              "DT","edgeR","org.Hs.eg.db","Glimma","magick","ComplexHeatmap",
              "plyr", "EnhancedVolcano","ggsci","gridExtra","ggpubr","rrvgo",
              "RColorBrewer","simplifyEnrichment","scatterplot3d","shiny","shinyjs","shinybusy",
              "shinycssloaders","viridis","shinymanager","shinythemes","hrbrthemes","InteractiveComplexHeatmap")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      BiocManager::install(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
