suppressPackageStartupMessages({
  library(biomaRt)
  library(BASiCS)
  library(scDD)
  library(BPSC)
  library(ggrepel)
  library(apeglm)
  library(heatmap3)
  library(ggplot2)
  library(ggfortify)
  library(stringr)
  library(RColorBrewer)
  library(MKmisc)
  library(DESeq2)
  library(Rtsne)
  library(MAST)
  library(reticulate)
  library(scde)
  library(edgeR)
  library(GGally)
  library(GSEABase)
  library(limma)
  library(reshape2)
  library(data.table)
  library(knitr)
  library(stringr)
  library(NMF)
  library(rsvd)
  library(RColorBrewer)
  library(MAST)
  library(pcaMethods)
  library(segmented)
  library(robust)
  library(MASS)
  library(umap)
  library(QoRTs)
  library(tximport)
  library(Seurat)
  require(ReactomePA)
  require(clusterProfiler)
  require(org.Sc.sgd.db)
  require(meshes)
  library(msigdbr)
  library(Seurat)
  library(doParallel)
  library(grid)
  library(gridExtra)
  library(pathview)
  library(igraph)
  library(ggraph)
  library(qvalue)
  library(ggrepel)
  library(DRIMSeq)
  library(ggbeeswarm)
  library(stageR)
  library(DEXSeq)
  library(trackViewer)
  library(GenomicFeatures)
})

registerDoParallel(cores=1)
#source('/Volumes/GoogleDrive/My Drive/Yeast/Analysis/FINAL_DATA/scripts/QCFilteringFunctions.R')
#source('/Volumes/GoogleDrive/My Drive/Yeast/Analysis/FINAL_DATA/scripts/Plots.R')

DOT_COLOR <- c(
               "#293acc",
               
               "#B20000",

               
               "#008331",
               
               
               "#000000"
               )
names(DOT_COLOR ) <- c( 
                       'mouseEgg', 
                        
                       'mouseZygote',

                       'ratEgg',
                       
                       'ratZygote')
CONDITIONS_NAMES <- c('Mouse Egg', 
                      'Mouse Zygote', 
                      'Rat Egg', 
                      'Rat Zygote'
                      )
names(CONDITIONS_NAMES) <- c( 
  'mouseEgg', 
  
  'mouseZygote',
  
  'ratEgg',
  
  'ratZygote')

use_python("/Users/ysu13/opt/anaconda3/bin/python")
FILES <- 'D:/gDrive/Colon_cancer/Milciclib/'
setwd('/Volumes/GoogleDrive/My Drive/mouse_rat/mouse_rat_final/')

# Wiki
WPGENE_Mm <- read.gmt(paste(FILES, "./wikipathways-20211110-gmt-Mus_musculus.gmt", sep = "")) %>% tidyr::separate( c("name","version","wpid","org"), "%")
WPGENE_Mm <- WPGENE[,c("wpid", "gene", "name")]


# Reactome
# delegated to ReactomPA

#MSigDb
MSDGENE <- as.data.frame(msigdbr(species = "Mus musculus", category = "C2")[,c("gs_id", "entrez_gene","gs_name")])
MSDGENE <- MSDGENE[!grepl("CANCER", MSDGENE[,3]), ]

#Biocyc

GSEnv2df <- function(env0){
  df <- do.call(rbind, lapply(names(env0), FUN = function(x){y <- env0[[x]]; cbind(rep(x, length(y$genes)),  y$genes, rep(y$name, length(y$genes)))}))
  row.names(df) <- 1:nrow(df)
  colnames(df) <- c("gs.id", "gene", "name" )
  df[is.na(df[,3]),3] <- df[is.na(df[,3]),1]
  df <- df[df[,2] != "YPR144W",]
  df[, 2] <- sapply(df[, 2], function(x){
    get(x, org.Sc.sgdENTREZID)
  })
  return(as.data.frame(df))
}
BIOCYC <- GSEnv2df(iso_v_aa$biocyc.env)

########################################################## INSTALLATIONS #######################################################

BiocManager::install(c("biomaRt",
"BASiCS",
"scDD",
"BPSC",
"heatmap3",
"pheatmap",
"ggplot2",
"ggfortify",
"stringr",
"RColorBrewer",
"MKmisc",
"DESeq2",
"Rtsne",
"MAST",
"reticulate",
"edgeR",
"GGally",
"GSEABase",
"limma",
"reshape2",
"data.table",
"knitr",
"NMF",
"rsvd",
"RColorBrewer",
"pcaMethods",
"segmented",
"robust",
"MASS",
"umap",
"QoRTs",
"tximport",
"Seurat",
"ReactomePA",
"clusterProfiler",
"org.Sc.sgd.db",
"meshes",
"msigdbr",
"doParallel",
"grid",
"gridExtra",
"universalmotif",
"memes"))

install.packages("devtools")
install_github("nghiavtr/BPSC")
install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",repos=NULL,type="source")
install.packages('Cairo')
devtools::install_version('flexmix', '2.3-13')
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)


