filterSegments <- function(segments, min.seg.width) {
if (is.null(segments)) {
return(NULL)
}
if (min.seg.width<=0) {
return(segments)
}
replace.index <- which(width(segments) < min.seg.width)
repl.segments <- segments[replace.index]
keep.index <- which(width(segments) >= min.seg.width)
keep.segments <- segments[keep.index]
nearest.index <- nearest(repl.segments, keep.segments)
na.mask <- is.na(nearest.index)
nearest.index <- nearest.index[!na.mask]
replace.index <- replace.index[!na.mask]
if (length(nearest.index)>0) {
nearest.keep.segments <- keep.segments[nearest.index]
segments$state[replace.index] <- nearest.keep.segments$state
segments$mstate[replace.index] <- nearest.keep.segments$mstate
segments$pstate[replace.index] <- nearest.keep.segments$pstate
}
segments.df <- as.data.frame(segments)
segments.df <- collapseBins(segments.df, column2collapseBy='state', columns2drop=c('width'))
segments.filtered <- as(segments.df, 'GRanges')
seqlevels(segments.filtered) <- seqlevels(segments) # correct order after as()
seqlengths(segments.filtered) <- seqlengths(segments)
return(segments.filtered)
}