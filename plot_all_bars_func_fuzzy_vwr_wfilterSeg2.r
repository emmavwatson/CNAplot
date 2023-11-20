#requires AneuFinder
#requires ggplot2
#requires subcloneCorrection_after_filterSeg2
#requires transCoord
#requires filterSegments


#chromColors <- c("brown1","cyan4", "darkgoldenrod2", "darkorange1","cadetblue2","brown1","cyan4","darkgoldenrod2","darkorange1","cadetblue2","brown1" ,"cyan4","darkgoldenrod2","darkorange1","cadetblue2","brown1","cyan4" ,"darkgoldenrod2","darkorange1","cadetblue2","brown1","cyan4","darkgoldenrod2")

#plot_all_bars_func_fuzzy_vwr_wfilterSeg2(returned_062119_plate2_dips_filterSeg4_fwa, 2, chromColors, 0, 4, files_062119_plate2_dips_filterSeg, 10000000, 0.7, "becool", 0.96, 3.05, 0.05)




plot_all_bars_func_fuzzy_vwr_wfilterSeg2 <- function(x,w,u,t,m,files,msw,maxvar,id,low,high,maxdif) {
    for(j in 1:length(x)) {
    rm(modelgr)
    rm(files_i)
    rm(dfplot.seg)
    rm(dfplot.segX)
    rm(binsmodelgrTC)
    rm(binsmodelgr)
    rm(bp.coords)
    rm(dfplot.segXnew)
    rm(plot1)
    rm(counts)
    modelgr <- x[[j]]
    files_i <- basename(files[j])
    binsmodelgr <- modelgr$bins
    bp.coords <- modelgr$breakpoints
    modelgr$segments <- filterSegments(modelgr$segments, msw)
    dfplot.seg <- as.data.frame(transCoord(modelgr$segments))
    dfplot.seg <- dfplot.seg[!(dfplot.seg$copy.number==0),]
    dfplot.seg$copy.number <- substr(dfplot.seg$state, 1, 1)
    dfplot.seg$copy.number <- as.numeric(dfplot.seg$copy.number)
    dfplot.seg$counts.CNV <- modelgr$distributions[as.character(dfplot.seg$state),'mu']
    dfplot.seg$width <- dfplot.seg$end.genome - dfplot.seg$start.genome
    binsmodelgrTC <- transCoord(binsmodelgr)
    binsmodelgrTC <- as.data.frame(binsmodelgrTC)
    counts <- as.data.frame(modelgr$bincounts)
    binsmodelgrTC <- cbind(binsmodelgrTC, counts)
    colnames(binsmodelgrTC) <- make.unique(names(binsmodelgrTC))   
    dfplot.segX <- subcloneCorrection_after_filterSeg2(dfplot.seg, w, binsmodelgrTC,maxdif)  
    colnames(dfplot.segX) <- make.unique(names(dfplot.segX))
    rm(dataset1)
    rm(temp_dataset1)
    rm(chr)
    chr <- split(dfplot.segX, dfplot.segX$seqnames)
        for (i in 1:length(chr)) {
                if (!exists("dataset1")){
                    dataset1 <- data.frame(chr[[i]])
                if (nrow(dataset1) <3) {next}
                if (nrow(dataset1) >10) {
                    dataset1$copy.numberNew <- mean(dataset1$copy.numberNew) 
                    next}
                if(var(dataset1$copy.numberNew) > maxvar) {next}
                    dataset1$copy.numberNew <- mean(dataset1$copy.numberNew)
            }
    
            if (exists("dataset1")){
                temp_dataset1 <- data.frame(chr[[i]])
                if (nrow(temp_dataset1) <3) {
                    dataset1 <-rbind(dataset1, temp_dataset1)
                    next}
                if (nrow(temp_dataset1) >10) {
                    temp_dataset1$copy.numberNew <- mean(temp_dataset1$copy.numberNew)
                    dataset1 <-rbind(dataset1, temp_dataset1)
                    next}
                if(var(temp_dataset1$copy.numberNew) > maxvar) {
                    dataset1 <-rbind(dataset1, temp_dataset1)
                    next}
                temp_dataset1$copy.numberNew <- mean(temp_dataset1$copy.numberNew)
                dataset1 <-rbind(dataset1, temp_dataset1)
                rm(temp_dataset1)
            }
        }
    dfplot.segXnew <- dataset1
    plot1 <- ggplot(data = binsmodelgrTC, aes(x = (start.genome+end.genome)/2, y = counts/(sum(dfplot.segXnew$width[dfplot.segXnew$copy.number==w]*dfplot.segXnew$counts.CNV[dfplot.segXnew$copy.number==w])/sum(as.numeric(dfplot.segXnew$width[dfplot.segXnew$copy.number==w])))*w, colour = seqnames)) + geom_point(size = 1.2, alpha= 0.1) + theme_classic() + scale_color_manual(values = u) + geom_segment(data=dfplot.segXnew, mapping=aes_string(x='start.genome',y='copy.numberNew',xend='end.genome',yend='copy.numberNew'), size=4) + scale_x_continuous(breaks=seqlengths(modelgr$bins)/2+cum.seqlengths.0[as.character(seqlevels(modelgr$bins))], labels=seqlevels(modelgr$bins)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=2) + ylim(t, m) + ylab("copy number") + xlab("chromosome") + theme(legend.position="none") + theme(plot.margin=grid::unit(c(0.5,0.2,0.5,0.2), "cm")) + theme(axis.text.y=element_text(size=16), axis.text.x=element_text(size=13), axis.title = element_text(size=13))
    ggsave(plot1, file=paste0("model_fuzzy_vwr_",files_i,id,"_barplot",".pdf"), height = 4.15, width = 8)
    write.csv(dfplot.segXnew, file = paste("segments", files_i, id, ".csv",sep=""), row.names = FALSE, col.names = FALSE)
    plot2 <- ggplot(data = binsmodelgrTC, aes(x = (start.genome+end.genome)/2, y = counts/(sum(dfplot.segXnew$width[dfplot.segXnew$copy.number==w]*dfplot.segXnew$counts.CNV[dfplot.segXnew$copy.number==w])/sum(as.numeric(dfplot.segXnew$width[dfplot.segXnew$copy.number==w])))*w)) + geom_point(size = 1.2, alpha= 0) + theme_classic() + geom_segment(data=dfplot.segXnew, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4, color='copy.numberNew'), size=6.5) + geom_segment(data=whiteout, mapping=aes_string(x='start.genome',y=4, xend='end.genome',yend=4), color="white", size=6.5) + scale_x_continuous(breaks=seqlengths(modelgr$bins)/2+cum.seqlengths.0[as.character(seqlevels(modelgr$bins))], labels=seqlevels(modelgr$bins)) + geom_vline(aes_string(xintercept='x'), data=df.chroms, col='black', linetype=1) + ylim(4,4)  + theme(legend.position="none") + theme(plot.margin=grid::unit(c(5,0,5,0), "cm")) + ylab("") + xlab("") + scale_colour_gradient2(limits = c(low, high), low="blue", mid="gray90", high="red", na.value="firebrick", midpoint=w) + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.line.y=element_blank()) 
    ggsave(plot2, file=paste0("model_bar_vwr_",files_i,id,"_barplot",".pdf"), height = 4.15, width = 8)
    } 
    }