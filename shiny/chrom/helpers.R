annotate_calls <- function(calls,intens=g_intens,intens_rev=g_intens_rev){
    iG <- intens[calls[,trace_peak]][,G]
    iA <- intens[calls[,trace_peak]][,A]
    iT <- intens[calls[,trace_peak]][,T]
    iC <- intens[calls[,trace_peak]][,C]

    if(!is.null(intens_rev)){
        iG_rev <- intens_rev[calls[,trace_peak]][,G]
        iA_rev <- intens_rev[calls[,trace_peak]][,A]
        iT_rev <- intens_rev[calls[,trace_peak]][,T]
        iC_rev <- intens_rev[calls[,trace_peak]][,C]
    }else{
        iG_rev <- iA_rev <- iT_rev <- iC_rev <- 0
    }

    calls[,c("iG_fwd","iA_fwd","iT_fwd","iC_fwd","iG_rev","iA_rev","iT_rev","iC_rev"):= list(iG,iA,iT,iC,iG_rev,iA_rev,iT_rev,iC_rev)]

    #contains codons table
    load("../../data/codons.rdata")
    calls       <-  merge(x = calls, y = codons[,list(gen_coord,codon,ord_in_cod,coding_seq,AA)], by = "gen_coord", all.x = TRUE)
    return(calls)
}

complement <- function(base){
    return (chartr("ATGC","TACG",base))
}
