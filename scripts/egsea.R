if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("PADOG", "GSVA", "AnnotationDbi", "topGO",
                       "pathview", "gage", "globaltest", "limma", "edgeR", "safe",
                       "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db"))




gs.annots = buildIdx(entrezIDs=row.names(test_limma$E), species="human", 
                     msigdb.gsets=c("h"), go.part = TRUE)


gsets <- buildIdx(
  row.names(test_limma$E),
  species = "human", msigdb.gsets = c("h"), kegg.exclude = c("Metabolism"))


test_limma <- egsea_limma_obj(counts(salmon_res5$dds), colData(salmon_res5$dds), control='control')
msiggs <- buildMSigDBIdx(
  row.names(test_limma$E),
  species = "Homo sapiens",
  geneSets = "h",
  go.part = FALSE,
  min.size = 10
)
contrast.matrix = makeContrasts('grouptreated-groupcontrol',
                                levels=test_limma$design)
gsa = egsea(voom.results = test_limma, contrasts = contrast.matrix, gs.annots = msiggs,
            symbolsMap = test_limma$genes, baseGSEAs =  egsea.base(), report.dir = "./egsea_res",
            sort.by = "p.adj", num.threads = 4, report = T, interactive = T, 
            keep.base=T, keep.set.scores = T)




