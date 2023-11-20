remove_nonconsecutive_outliers2 <- function(x,y) {  
    bins_outrm <- x
    bins_outrm$outliers <- 0
    shiftn <- function(x, n){
        c(x[-(seq(n))], rep(NA, n)) }
    for (j in 1:22) {
        outliers <- boxplot(subset(bins_outrm, seqnames==j, select =c(counts)), plot=FALSE, range=y)$out
        bins_outrm$outliers <- ifelse(bins_outrm$seqnames==j & (bins_outrm$counts %in% outliers), "outlier", bins_outrm$outliers)
    }
    outliers <- boxplot(subset(bins_outrm, seqnames=="X", select =c(counts)), plot=FALSE, range=y)$out
    bins_outrm$outliers <- ifelse(bins_outrm$seqnames=="X" & (bins_outrm$counts %in% outliers), "outlier", bins_outrm$outliers)
    bins_outrm$pre <- bins_outrm$outliers
    bins_outrm$pre <- shiftn(bins_outrm$pre, 1)
    bins_outrm$pre_counts <- shiftn(bins_outrm$counts, 1)
    bins_outrm$lonely_outlier <- ifelse(bins_outrm$outliers =="outlier" & bins_outrm$pre!="outlier", 1, ifelse(bins_outrm$outliers =="outlier" & bins_outrm$pre=="outlier" & ((bins_outrm$counts/bins_outrm$pre_counts) < 0.7 | (bins_outrm$counts/bins_outrm$pre_counts) > 1.3),1,0))
    bins_outrm <- subset(bins_outrm, lonely_outlier=="0")
    bins_outrm <- bins_outrm[,c(1:6)]
    return(bins_outrm) }