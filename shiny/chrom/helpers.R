annotate_calls <- function(calls,intens){
    iG <- intens[calls[,trace_peak]][,G]
    iA <- intens[calls[,trace_peak]][,A]
    iT <- intens[calls[,trace_peak]][,T]
    iC <- intens[calls[,trace_peak]][,C]
    calls[,c("G","A","T","C"):= list(iG,iA,iT,iC)]
    
    cod         <- load("../../data/codons.rdata")
    
    calls       <-  merge(x = calls, y = cod[,list(gen_coord,cod,ord_in_cod)], by = "gen_coord", all.x = TRUE)
    return(calls)
}

complement <- function(base){
    return (chartr("ATGC","TACG",base))
}