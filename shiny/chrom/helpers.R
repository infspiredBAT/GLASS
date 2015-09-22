annotate_calls <- function(calls,intens,intens_rev){

    #contains codons table
    load("../../data/codons.rdata")
    calls <- merge(x = calls, y = codons[,list(gen_coord,codon,ord_in_cod,coding_seq,AA_ref=AA,AA_mod=AA)], by = "gen_coord", all.x = TRUE)
    #reorder columns so that id is first (so that the checkbox from the shiny data table selects the correct value and the delete button knows what to delete)
    setcolorder(calls,c("id",colnames(calls)[-2]))

    #calculate the noise levels
    calls[,noise_abs_fwd:=noise(iA_fwd,iC_fwd,iG_fwd,iT_fwd,TRUE),by=1:nrow(calls)]
    calls[,noise_rel_fwd:=noise(iA_fwd,iC_fwd,iG_fwd,iT_fwd,FALSE),by=1:nrow(calls)]
    if("call_rev" %in% colnames(calls)){
        calls[,noise_abs_rev:=noise(iA_rev,iC_rev,iG_rev,iT_rev,TRUE),by=1:nrow(calls)]
        calls[,noise_rel_rev:=noise(iA_rev,iC_rev,iG_rev,iT_rev,FALSE),by=1:nrow(calls)]
    }

    #precalculate neighbourhood  (absolute,relative)x(forward,reverse)
    nbrhd_a_f <- rollmean(calls$noise_abs_fwd,k=7)
    nbrhd_r_f <- rollmean(calls$noise_rel_fwd,k=7)
    if("call_rev" %in% colnames(calls)){
        nbrhd_a_r <- rollmean(calls$noise_abs_rev,k=7)
        nbrhd_r_r <- rollmean(calls$noise_rel_rev,k=7)
    }

    calls[,rm7noise_abs_fwd := c(rep(nbrhd_a_f[1],3),nbrhd_a_f,rep(nbrhd_a_f[length(nbrhd_a_f)],3))]
    calls[,rm7noise_rel_fwd := c(rep(nbrhd_r_f[1],3),nbrhd_r_f,rep(nbrhd_r_f[length(nbrhd_r_f)],3))]
    if("call_rev" %in% colnames(calls)){
        calls[,rm7noise_abs_rev := c(rep(nbrhd_a_r[1],3),nbrhd_a_r,rep(nbrhd_a_r[length(nbrhd_a_r)],3))]
        calls[,rm7noise_rel_rev := c(rep(nbrhd_r_r[1],3),nbrhd_r_r,rep(nbrhd_r_r[length(nbrhd_r_r)],3))]
    }

    #ref peak = pseudo trace = sum of intensities
    calls[,ref_peak_abs_fwd:=sum(iA_fwd,iC_fwd,iG_fwd,iT_fwd),by=1:nrow(calls)]
    if("call_rev" %in% colnames(calls)) calls[,ref_peak_abs_rev:=sum(iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]

    #add information about the first and second highest peak
    calls[,c("sample_peak_base_fwd","sample_peak_abs_fwd") := i_wo_p(1,iA_fwd,iC_fwd,iG_fwd,iT_fwd),by=1:nrow(calls)]
    calls[,sample_peak_pct_fwd := ((100/ref_peak_abs_fwd)*sample_peak_abs_fwd)]
    calls[,c("mut_peak_base_fwd","mut_peak_abs_fwd") := i_wo_p(2,iA_fwd,iC_fwd,iG_fwd,iT_fwd),by=1:nrow(calls)]
    calls[,mut_peak_pct_fwd := ((100/ref_peak_abs_fwd)*mut_peak_abs_fwd)]
    calls[,mut_s2n_abs_fwd:=mut_peak_abs_fwd/noise_abs_fwd]
    calls[,mut_call_fwd:=sample_peak_base_fwd]
    if("call_rev" %in% colnames(calls)){
        calls[,c("sample_peak_base_rev","sample_peak_abs_rev") := i_wo_p(1,iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]
        calls[,sample_peak_pct_rev := ((100/ref_peak_abs_rev)*sample_peak_abs_rev)]
        calls[,c("mut_peak_base_rev","mut_peak_abs_rev") := i_wo_p(2,iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]
        calls[,mut_peak_pct_rev := ((100/ref_peak_abs_rev)*mut_peak_abs_rev)]
        calls[,mut_s2n_abs_rev:=mut_peak_abs_rev/noise_abs_rev]
        calls[,mut_call_rev:=sample_peak_base_rev]
    }
    calls[,set_by_user:=FALSE]
    calls[,cons:=user_sample]
    return(calls)
}

call_variants <- function(calls, qual_thres, mut_min, s2n_min){
    # reset all but set_by_user
    calls[set_by_user == FALSE, user_sample := cons]
    calls[set_by_user == FALSE, user_mut := cons]

    # calls[(rm7qual < qual_thres | quality < qual_thres) & set_by_user == FALSE, user_sample := "low qual"]

    # mut
    if("call_rev" %in% colnames(calls)) {
        calls[    mut_peak_pct_fwd >= mut_min
                & mut_s2n_abs_fwd >= s2n_min
                # & abs(mut_peak_pct_fwd - mut_peak_pct_rev) < mut_min
                # & mut_peak_base_fwd == mut_peak_base_rev
                & !set_by_user,
                c("mut_call_fwd", "user_mut") := list(mut_peak_base_fwd, mut_peak_base_fwd)
                # user_mut := mut_call_fwd
                ]
        calls[    mut_peak_pct_rev >= mut_min
                & mut_s2n_abs_rev >= s2n_min
                & mut_peak_pct_rev > mut_peak_pct_fwd
                & !set_by_user,
                c("mut_call_rev", "user_mut") := list(mut_peak_base_rev, mut_peak_base_rev)
                ]
    } else {
        calls[    mut_peak_pct_fwd >= mut_min
                & mut_s2n_abs_fwd >= s2n_min
                & !set_by_user,
                #user_sample := paste(user_sample,tolower(mut_peak_base_fwd),sep="")]
                user_mut := mut_call_fwd
            ]
    }
    return(calls)
}

complement <- function(base){
    return (chartr("ATGC","TACG",base))
}
retranslate <- function(calls){
    return(calls[coding_seq>0,AA_mod := rep(translate(calls[coding_seq>0,user_sample],frame=calls[coding_seq>0,ord_in_cod][1] -1),each=3)])
}

#background noise absolute or relative to reference peak
noise <- function(a,b,c,d,abs=FALSE){
    vec <- sort(c(a,b,c,d))
    if(a==0 & b==0 & c==0 & d==0) return(0)
    if(abs) return(mean(vec[1:3]))
    else    return(mean(vec[1:3]/sum(vec[4])))
}

#retrieve the name, intens value, and with/without trace pos of the p'th peak as this may differ from the "call" which may bey ambiguous iupac (S,W,R...)
i_wo_p <- function(p,iA,iC,iG,iT){
    mut_peak <- (sort(c(iA,iC,iG,iT),decreasing = TRUE)[p])
         if (mut_peak == 0 ) return(list(" ",mut_peak))
    else if (mut_peak == iA) return(list("A",mut_peak))
    else if (mut_peak == iC) return(list("C",mut_peak))
    else if (mut_peak == iG) return(list("G",mut_peak))
    else if (mut_peak == iT) return(list("T",mut_peak))
}
i_w_p <- function(p,iA,iC,iG,iT,pA,pC,pG,pT){
    mut_peak <- (sort(c(iA,iC,iG,iT),decreasing = TRUE)[p])
         if (mut_peak == 0 ) return(list(" ",mut_peak,pA))
    else if (mut_peak == iA) return(list("A",mut_peak,pA))
    else if (mut_peak == iC) return(list("C",mut_peak,pC))
    else if (mut_peak == iG) return(list("G",mut_peak,pG))
    else if (mut_peak == iT) return(list("T",mut_peak,pT))
}
