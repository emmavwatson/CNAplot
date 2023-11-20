#files_scDNA_test <- list.files("/path/to/seqfile/BAMs/folder", full.names = TRUE)

#Returned_models_files <- lapply(files_scDNA_test, function(x) {
#bedfile <- file.path(x)
#model <- final_plot_func3_modelreturn_cluster_vwr_dips_outSalvage5(bedfile, 3e6, 0.7, blacklist_strict3, BSgenome.Hsapiens.1000genomes.hs37d5, "hg19", vwr)
#return(model)})

#models can be plotted with feature extraction using: plot_all_bars_func_fuzzy_vwr_wfilterSeg2.r




final_plot_func3_modelreturn_cluster_vwr_dips_outSalvage5 <- function(x,y,t,f,BS,assembly, vwr) {
    bins <- binReads(x, assembly=assembly, binsizes=y, chromosomes=c(1:22, 'X'), blacklist=f, variable.width.reference = vwr)
    binsGC <- correctGC(list(bins[[1]]), GC.BSgenome = BS)
    grbinsGC <-unlist(binsGC[[1]])
    df <- data.frame(seqnames=seqnames(grbinsGC),start=start(grbinsGC)-1,end=end(grbinsGC),names=c(rep(".", length(grbinsGC))),counts=elementMetadata(grbinsGC)$counts,strands=strand(grbinsGC))
    df <- subset(df, counts > 1)
    dfbins_rmout <- remove_nonconsecutive_outliers2(df, t)
    seginfo <- as.data.frame(seqlengths(grbinsGC))
    seginfo$seqnames <- rownames(seginfo)
    seqinfoG23_C1_25 <- Seqinfo(seqnames=seginfo[,2], seqlengths=seginfo[,1], isCircular = NA, genome = NA)
    gr <- makeGRangesFromDataFrame(dfbins_rmout,keep.extra.columns=TRUE,ignore.strand=FALSE, seqinfo=seqinfoG23_C1_25, seqnames.field=c("seqnames"), start.field="start", end.field=c("end"), strand.field="strand", starts.in.df.are.0based=FALSE)
    modelgr <- findCNVs(gr, method = "HMM", most.frequent.state = '2-somy', max.iter = 80)
    return(modelgr) }
