library("org.Hs.eg.db")
library(GOstats)
library(Category)
library(qvalue)
analyze <- function(ens, df){
  clean_and_convert <- function(ensembl_ids) {
    cleaned <- sub("\\..*", "", ensembl_ids)
    cleaned <- cleaned[!is.na(cleaned)]
    cleaned <- cleaned[nchar(cleaned) > 0]
    suppressMessages(
      mapIds(org.Hs.eg.db, 
             keys = cleaned,
             column = "ENTREZID",
             keytype = "ENSEMBL",
             multiVals = "first")
    )
  }
  sel_ensembl <- unique(na.omit(ens))
  if(length(sel_ensembl) == 0) stop("No valid selected ENSEMBL IDs")
  sel_entrez <- clean_and_convert(sel_ensembl)
  sel_entrez <- na.omit(sel_entrez)
  if(length(sel_entrez) == 0) stop("No valid Entrez IDs in selection")
  all_ensembl <- unique(na.omit(df$name))
  if(length(all_ensembl) == 0) stop("No valid background ENSEMBL IDs")
  all_entrez <- clean_and_convert(all_ensembl)
  all_entrez <- na.omit(all_entrez)
  if(length(all_entrez) == 0) stop("No valid Entrez IDs in background")
  tryCatch({
    params <- new("GOHyperGParams",
                  geneIds = sel_entrez,
                  universeGeneIds = all_entrez,
                  ontology = "BP",
                  pvalueCutoff = 0.05,
                  conditional = TRUE,
                  testDirection = "over",
                  annotation = "org.Hs.eg.db")
    
    over <- hyperGTest(params)
    ov <- summary(over)
    ov <- ov[ov[,6] <= 500 & ov[,6] >= 3, ]
    ov <- cbind(ov, 
                idx = sapply(ov$GOBPID, function(go_id) {
                  genes <- get(go_id, org.Hs.egGO2ALLEGS)
                  paste(intersect(genes, sel_entrez), collapse = ", ")
                }),
                symbol = sapply(ov$GOBPID, function(go_id) {
                  genes <- get(go_id, org.Hs.egGO2ALLEGS)
                  symbols <- mapIds(org.Hs.eg.db, 
                                    keys = intersect(genes, sel_entrez),
                                    column = "SYMBOL",
                                    keytype = "ENTREZID")
                  paste(symbols, collapse = ", ")
                }))
    ov$p_adjusted <- p.adjust(ov$Pvalue, "BH")
    
    return(ov)
  }, error = function(e) {
    message("Error in GO analysis: ", e$message)
    return(data.frame())
  })
}



