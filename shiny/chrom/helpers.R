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
    calls <- merge(x = calls, y = codons[,list(gen_coord,codon,ord_in_cod,coding_seq,AA)], by = "gen_coord", all.x = TRUE)
    #reorder columns so that id is first (so that the checkbox from the shiny data table selects the correct value and the delete button knows what to delete)
    setcolorder(calls,c("id",colnames(calls)[-2]))

    #calculate the noise levels
    calls[,noise_abs_fwd:=noise(iG_fwd,iC_fwd,iA_fwd,iT_fwd,TRUE),by=1:nrow(calls)]
    calls[,noise_rel_fwd:=noise(iG_fwd,iC_fwd,iA_fwd,iT_fwd,FALSE),by=1:nrow(calls)]
    if(!is.null(intens_rev)){
        calls[,noise_abs_rev:=noise(iG_rev,iC_rev,iA_rev,iT_rev,TRUE),by=1:nrow(calls)]
        calls[,noise_rel_rev:=noise(iG_rev,iC_rev,iA_rev,iT_rev,FALSE),by=1:nrow(calls)]
    }

    #precalculate neighbourhood  (absolute,relative)x(forward,reverse)
    nbrhd_a_f <- rollmean(calls$noise_abs_fwd,k=7)
    nbrhd_r_f <- rollmean(calls$noise_rel_fwd,k=7)
    if(!is.null(intens_rev)){
        nbrhd_a_r <- rollmean(calls$noise_abs_rev,k=7)
        nbrhd_r_r <- rollmean(calls$noise_rel_rev,k=7)
    }

    calls[,rm7noise_abs_fwd := c(rep(nbrhd_a_f[1],3),nbrhd_a_f,rep(nbrhd_a_f[length(nbrhd_a_f)],3))]
    calls[,rm7noise_rel_fwd := c(rep(nbrhd_r_f[1],3),nbrhd_r_f,rep(nbrhd_r_f[length(nbrhd_r_f)],3))]
    if(!is.null(intens_rev)){
        calls[,rm7noise_abs_rev := c(rep(nbrhd_a_r[1],3),nbrhd_a_r,rep(nbrhd_a_r[length(nbrhd_a_r)],3))]
        calls[,rm7noise_rel_rev := c(rep(nbrhd_r_r[1],3),nbrhd_r_r,rep(nbrhd_r_r[length(nbrhd_r_r)],3))]
    }

    calls[,ref_peak_abs_fwd:=sum(iG_fwd,iA_fwd,iT_fwd,iC_fwd),by=1:nrow(calls)]
    if(!is.null(intens_rev)) calls[,ref_peak_abs_rev:=sum(iG_rev,iA_rev,iT_rev,iC_rev),by=1:nrow(calls)]

    #add information about the second highest peak
    calls[,c("mut_peak_call_fwd","mut_peak_abs_fwd"):=second(iA_fwd,iC_fwd,iG_fwd,iT_fwd),by=1:nrow(calls)]
    calls[,mut_s2n_abs_fwd:=mut_peak_abs_fwd/noise_abs_fwd]
    calls[,mut_peak_pct_fwd:=((100/ref_peak_abs_fwd)*mut_peak_abs_fwd)]
    if(!is.null(intens_rev)){
        calls[,c("mut_peak_call_rev","mut_peak_abs_rev"):=second(iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]
        calls[,mut_s2n_abs_rev:=mut_peak_abs_rev/noise_abs_rev]
        calls[,mut_peak_pct_rev:=((100/ref_peak_abs_rev)*mut_peak_abs_rev)]
    }
    calls[,reset_by_user:=FALSE]
    calls[,cons:=user_mod]
    return(calls)
}

call_variants <- function(calls, qual_thres, mut_min, s2n_min){

    # reset all but reset_by_user
    calls[reset_by_user == FALSE, user_mod := cons]

    # low qual
    calls[(rm7qual < qual_thres | quality < qual_thres) & reset_by_user == FALSE, user_mod := "low qual"]

    # mut peak
    if("call_rev" %in% colnames(calls)) {
        calls[(mut_peak_pct_fwd > mut_min & mut_peak_pct_rev > mut_min & mut_s2n_abs_fwd > s2n_min & mut_s2n_abs_rev > s2n_min & mut_peak_call_fwd == mut_peak_call_rev & user_mod != "low qual"), user_mod := paste(user_mod,mut_peak_call_fwd,sep="")]
    } else {
        calls[(mut_peak_pct_fwd > mut_min                              & mut_s2n_abs_fwd > s2n_min                                                                      & user_mod != "low qual"), user_mod := paste(user_mod,mut_peak_call_fwd,sep="")]
    }

    return(calls)
}

complement <- function(base){
    return (chartr("ATGC","TACG",base))
}

#background noise absolute or relative to reference peak
noise <- function(a,b,c,d,abs=FALSE){
    vec <- sort(c(a,b,c,d))
    if(a==0 & b==0 & c==0 & d==0) return(0)
    if(abs) return(mean(vec[1:3]))
    else    return(mean(vec[1:3]/sum(vec[4])))
}

#retrieve the name and value of the second highest peak
second <- function(iA,iC,iG,iT){
    mut_peak <- (sort(c(iA,iC,iG,iT),decreasing = TRUE)[2])
         if (mut_peak == iA) return(list("A",mut_peak))
    else if (mut_peak == iC) return(list("C",mut_peak))
    else if (mut_peak == iG) return(list("G",mut_peak))
    else if (mut_peak == iT) return(list("T",mut_peak))
}
