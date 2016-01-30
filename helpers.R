annotate_calls <- function(calls,intens,intens_rev,glassed_cod){

    #contains codons table
    load(glassed_cod)
    calls <- merge(x = calls, y = cod_table[,list(gen_coord,codon,ord_in_cod,coding_seq,aa_ref=AA)], by = "gen_coord", all.x = TRUE)
    calls[aa_ref != "",aa_ref:=aaa(toupper(aa_ref))]
    calls[,c("aa_sample","aa_mut"):=aa_ref]
    cod_table <<- cod_table
    #reorder columns so that id is first (so that the checkbox from the shiny data table selects the correct value and the delete button knows what to delete)
    setcolorder(calls,c("id",colnames(calls)[-2]))

    #calculate the noise levels
    calls[,noise_abs_fwd:=noise(iA_fwd,iC_fwd,iG_fwd,iT_fwd,TRUE),by=1:nrow(calls)]
    # calls[,noise_rel_fwd:=noise(iA_fwd,iC_fwd,iG_fwd,iT_fwd,FALSE),by=1:nrow(calls)]
    calls[,roll_noise_abs_mean_ratio_fwd := rollapply(calls$noise_abs_fwd, 11, function(x) (mean(x[1:5])+0)/(mean(x[7:11])+0), fill = "1")]
    # calls[,roll_noise_rel_mean_ratio_fwd := rollapply(calls$noise_rel_fwd, 11, function(x) (mean(x[1:5])+1)/(mean(x[7:11])+1), fill = "1")]
    if("call_rev" %in% colnames(calls)){
        calls[,noise_abs_rev:=noise(iA_rev,iC_rev,iG_rev,iT_rev,TRUE),by=1:nrow(calls)]
        # calls[,noise_rel_rev:=noise(iA_rev,iC_rev,iG_rev,iT_rev,FALSE),by=1:nrow(calls)]
        calls[,roll_noise_abs_mean_ratio_rev := rollapply(calls$noise_abs_rev, 11, function(x) (mean(x[1:5])+0)/(mean(x[7:11])+0), fill = "1")]
        # calls[,roll_noise_rel_mean_ratio_rev := rollapply(calls$noise_rel_rev, 11, function(x) (mean(x[1:5])+1)/(mean(x[7:11])+1), fill = "1")]
    }

#     #precalculate neighbourhood  (absolute,relative)x(forward,reverse)
#     nbrhd_a_f <- rollmean(calls$noise_abs_fwd,k=7)
#     nbrhd_r_f <- rollmean(calls$noise_rel_fwd,k=7)
#     if("call_rev" %in% colnames(calls)){
#         nbrhd_a_r <- rollmean(calls$noise_abs_rev,k=7)
#         nbrhd_r_r <- rollmean(calls$noise_rel_rev,k=7)
#     }
#
#     calls[,rm7noise_abs_fwd := c(rep(nbrhd_a_f[1],3),nbrhd_a_f,rep(nbrhd_a_f[length(nbrhd_a_f)],3))]
#     calls[,rm7noise_rel_fwd := c(rep(nbrhd_r_f[1],3),nbrhd_r_f,rep(nbrhd_r_f[length(nbrhd_r_f)],3))]
#     if("call_rev" %in% colnames(calls)){
#         calls[,rm7noise_abs_rev := c(rep(nbrhd_a_r[1],3),nbrhd_a_r,rep(nbrhd_a_r[length(nbrhd_a_r)],3))]
#         calls[,rm7noise_rel_rev := c(rep(nbrhd_r_r[1],3),nbrhd_r_r,rep(nbrhd_r_r[length(nbrhd_r_r)],3))]
#     }

    #ref peak = pseudo trace = sum of intensities
    calls[,ref_peak_abs_fwd:=sum(iA_fwd,iC_fwd,iG_fwd,iT_fwd),by=1:nrow(calls)]
    if("call_rev" %in% colnames(calls)) calls[,ref_peak_abs_rev:=sum(iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]

    #add information about the first and second highest peak
    calls[,c("sample_peak_base_fwd","sample_peak_abs_fwd") := i_wo_p(1,iA_fwd,iC_fwd,iG_fwd,iT_fwd),by=1:nrow(calls)]
    calls[,sample_peak_pct_fwd := ((100/ref_peak_abs_fwd)*sample_peak_abs_fwd)]
    calls[,c("mut_peak_base_fwd","mut_peak_abs_fwd") := i_wo_p(2,iA_fwd,iC_fwd,iG_fwd,iT_fwd),by=1:nrow(calls)]
    calls[,mut_peak_pct_fwd := ((100/ref_peak_abs_fwd)*mut_peak_abs_fwd)]
    calls[,mut_s2n_abs_fwd:=mut_peak_abs_fwd/noise_abs_fwd]
    calls[,mut_call_fwd:=call]
    if("call_rev" %in% colnames(calls)){
        calls[,c("sample_peak_base_rev","sample_peak_abs_rev") := i_wo_p(1,iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]
        calls[,sample_peak_pct_rev := ((100/ref_peak_abs_rev)*sample_peak_abs_rev)]
        calls[,c("mut_peak_base_rev","mut_peak_abs_rev") := i_wo_p(2,iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]
        calls[,mut_peak_pct_rev := ((100/ref_peak_abs_rev)*mut_peak_abs_rev)]
        calls[,mut_s2n_abs_rev:=mut_peak_abs_rev/noise_abs_rev]
        calls[,mut_call_rev:=call_rev]

        suppressWarnings(calls[,sample_peak_pct := max(c(sample_peak_pct_fwd,sample_peak_pct_rev),na.rm=TRUE),by=1:nrow(calls)])
        calls[,mut_peak_pct := mean(c(mut_peak_pct_fwd,mut_peak_pct_rev),na.rm=TRUE),by=1:nrow(calls)]
    } else {
        calls[,sample_peak_pct := sample_peak_pct_fwd]
        calls[,mut_peak_pct := mut_peak_pct_fwd]
    }
    calls[,set_by_user:=FALSE]
    # calls[,user_sample_orig:=user_sample]
    calls[set_by_user == FALSE, user_sample_orig := ambig_minus(user_sample,reference),by=1:nrow(calls[set_by_user==FALSE,])]
    return(calls)
}

get_noisy_neighbors <- function(calls){
    if("call_rev" %in% colnames(calls)){
        noisy_neighbors <- calls[trace_peak != "NA" & !is.na(gen_coord) & call != "-" & call_rev != "-"
                                 & roll_noise_abs_mean_ratio_fwd > 0 & roll_noise_abs_mean_ratio_fwd < 10 & roll_noise_abs_mean_ratio_rev > 0 & roll_noise_abs_mean_ratio_rev < 10
                                 & roll_noise_abs_mean_ratio_fwd / roll_noise_abs_mean_ratio_rev >= 5
                                ]
    } else {
        noisy_neighbors <- calls[trace_peak != "NA" & !is.na(gen_coord) & call != "-"
                                 & roll_noise_abs_mean_ratio_fwd > 0 & roll_noise_abs_mean_ratio_fwd < 10
                                 & roll_noise_abs_mean_ratio_fwd >= 2
                                ]
    }
    setkey(noisy_neighbors,id)
    return(noisy_neighbors)
}

splice_variants <- function(intrexdat){
    intrexdat$intrex$splicevar <- mapply(function(x,y) if (x == 'intron9' & abs(133 - y) <= 1) { '|beta variant' } else { '' }, intrexdat$intrex$attr, intrexdat$intrex$length)
    return(intrexdat)
}


include_locked_indels <- function(calls,vec,indels,fwd){
    vec <- copy(vec)
    get_del_positions <- function(code,pos,vec){
        coord1 <- as.numeric(gsub("c\\.(\\d*).*","\\1",code))
        coord2 <- suppressWarnings(as.numeric(gsub("c\\.\\d*_(\\d*).*","\\1",code)))
        if(is.na(coord2)) del_len <- 0
        else del_len <- coord2 - coord1
        if(vec[pos + 1] == "-") pos <- pos + 1
        return(0:del_len + pos)
    }

    get_ins_positions <- function(code,pos,vec){
        if(length(grep("nt",code)) > 0) ins_len <- as.numeric(gsub(".*nt(\\d*)","\\1",code))
        else ins_len <- nchar(gsub(".*[ins|dup](.*)","\\1",code))
        if(all(calls$reference[pos] == "-")) pos <- pos - 1

        return(1:ins_len + floor(pos))
    }

    dels <- lapply(names(indels)[grep("del",names(indels))],function(x) get_del_positions(x,indels[[x]],vec))
    if(length(dels) > 0){
        sec_dels <- dels[which(sapply(dels,function(x) all(vec[x] != "-")))]
        prim_dels <- setdiff(dels,sec_dels)
    } else{
        sec_dels <- list()
        prim_dels <- list()
    }


    ins <- lapply(names(indels)[grep("ins|dup",names(indels))],function(x) get_ins_positions(x,indels[[x]],vec))
    if(length(ins) > 0){
        sec_ins <- ins[which(sapply(ins,function(x) all(calls$reference[x] != "-")))]
        prim_ins <- setdiff(ins,sec_ins)
    } else{
        sec_ins <- list()
        prim_ins <- list()
    }

    move_vec <- numeric(length(vec))
    if(fwd){
        if(length(sec_dels) > 0) move_vec[sapply(sec_dels,function(x) max(x) + 1)] <- - sapply(sec_dels,length)
        if(length(prim_dels) > 0) move_vec[sapply(prim_dels,min)] <- sapply(prim_dels,length)
        if(length(sec_ins) > 0) move_vec[sapply(sec_ins,min)] <- sapply(sec_ins,length)
        if(length(prim_ins) > 0) move_vec[sapply(prim_ins,function(x) max(x) + 1)] <- -sapply(prim_ins,length)
        move_vec <- cumsum(move_vec)
    } else {
        if(length(sec_dels) > 0) move_vec[sapply(sec_dels,function(x) min(x) - 1)] <- sapply(sec_dels,length)
        if(length(prim_dels) > 0) move_vec[sapply(prim_dels,max)] <- -sapply(prim_dels,length)
        if(length(sec_ins) > 0)move_vec[sapply(sec_ins,max)] <- -sapply(sec_ins,length)
        if(length(prim_ins) > 0) move_vec[sapply(prim_ins,function(x) min(x) - 1)] <- sapply(prim_ins,length)
        move_vec <- rev(cumsum(rev(move_vec)))
    }

    new_vec <- rep("-",length(vec))
    vec[unlist(prim_dels)] <- calls$reference[unlist(prim_dels)]
    pos_vec <- seq_along(new_vec) + move_vec
    pos_vec[pos_vec < 1] <- 1
    new_vec[pos_vec] <- vec
    new_vec <- new_vec[1:length(vec)]
    new_vec[unlist(prim_dels)] <- "-"
    if(fwd) new_vec[unlist(sec_ins)] <- calls$mut_peak_base_fwd[unlist(sec_ins)]
    else new_vec[unlist(sec_ins)] <- calls$mut_peak_base_rev[unlist(sec_ins)]

#     if(length(het_dels) > 0){
#         vec <- vec[-het_dels]
#         if(fwd) {
#             vec <- c(vec,rep("-",length(het_dels)))
#         } else {
#             vec <- c(rep("-",length(het_dels)),vec)
#         }
#     }
#
#     prim_dels <- setdiff(dels,het_dels)
#     if(length(prim_dels) > 0){
#         vec[prim_dels] <- g_calls$reference[prim_dels]
#         new_vec <- rep("-",length(vec))
#         if(fwd) {
#             vec <- vec[1:min(length(vec),length(new_vec) - length(prim_dels))]
#             new_vec[setdiff(seq_along(new_vec),prim_dels)] <- vec
#         } else {
#             vec <- vec[1 + length(vec) -  min(length(vec),length(new_vec) - length(prim_dels)):1]
#             new_vec[setdiff(seq_along(new_vec),prim_dels)] <- vec
#         }
#         vec <- new_vec
#     }




    return(new_vec)
}

call_variants <- function(calls, qual_thres, mut_min, s2n_min,stored_het_indels){
    # reset all but set_by_user
    calls[set_by_user == FALSE, user_sample := user_sample_orig]
    calls[set_by_user == FALSE, user_mut    := user_sample_orig]
    calls[set_by_user == FALSE, mut_call_fwd := call]
    if(length(grep("del|ins|dup",names(stored_het_indels))) > 0){
        #g_indels_present <<- TRUE
        calls[, mut_call_fwd := include_locked_indels(calls,mut_call_fwd,stored_het_indels,fwd = T)]
    }else#{
    #    g_indels_present <<- FALSE
    #}
    # calls[set_by_user == FALSE, mut_call_fwd := ambig_minus(call,reference),by=1:nrow(calls[set_by_user==FALSE,])]
    # mut
    if("call_rev" %in% colnames(calls)) {
        # reset all but set_by_user
        calls[set_by_user == FALSE, mut_call_rev := call_rev]
        if(length(grep("del|ins|dup",names(stored_het_indels))) > 0){
            calls[, mut_call_rev := include_locked_indels(mut_call_rev,stored_het_indels,fwd = F)]
        }

        # calls[set_by_user == FALSE, mut_call_rev := ambig_minus(call_rev,reference),by=1:nrow(calls[set_by_user==FALSE,])]
        # initialising mut calls
        calls[
              mut_peak_pct_fwd >= mut_min
            & mut_s2n_abs_fwd >= s2n_min
            #& quality_fwd >= qual_thres
            , mut_call_fwd := mut_peak_base_fwd
            ]
        calls[set_by_user == FALSE, mut_call_fwd := ambig_minus(mut_call_fwd,reference),by=1:nrow(calls[set_by_user==FALSE,])]
        calls[
              mut_peak_pct_rev >= mut_min
            & mut_s2n_abs_rev >= s2n_min
            #& quality_rev >= qual_thres
            , mut_call_rev := mut_peak_base_rev
            ]
        calls[set_by_user == FALSE, mut_call_rev := ambig_minus(mut_call_rev,reference),by=1:nrow(calls[set_by_user==FALSE,])]
#         calls[
#               set_by_user == FALSE
#             #& mut_call_fwd != call
#             , c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
#             ]
        # setting user muts based on reference and quality
        calls[
              mut_call_fwd!=reference
            & mut_call_rev==reference
            & quality_fwd > qual_thres
            & set_by_user == FALSE
            , c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
            ]
        calls[
              mut_call_fwd==reference
            & mut_call_rev!=reference
            & quality_rev > qual_thres
            & set_by_user == FALSE
            , c("user_mut","mut_peak_pct") := list(mut_call_rev,mut_peak_pct_rev)
            ]
        calls[
              mut_call_fwd!=reference
            & mut_call_rev!=reference
            & set_by_user == FALSE
            # & mut_call_rev != call_rev
            , c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
            ]
        calls[
              mut_call_fwd!=reference
            & mut_call_rev!=reference
            & quality_rev > quality_fwd
            & set_by_user == FALSE
            # & mut_call_rev != call_rev
            , c("user_mut","mut_peak_pct") := list(mut_call_rev,mut_peak_pct_rev)
            ]
    } else {
        calls[
              mut_peak_pct_fwd >= mut_min
            & mut_s2n_abs_fwd >= s2n_min
            & mut_peak_base_fwd != reference
            , mut_call_fwd := mut_peak_base_fwd
            ]
        calls[set_by_user == FALSE, mut_call_fwd := ambig_minus(mut_call_fwd,reference),by=1:nrow(calls[set_by_user==FALSE,])]
        calls[
              set_by_user == FALSE
            & mut_call_fwd != call
            , c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
            ]
    }
    # masking low quality
    calls[quality < qual_thres & set_by_user == FALSE, c("user_sample","user_mut") := "N"]
    return(calls)
}

complement <- function(base){
    return (chartr("ATGCRYKMBVDH","TACGYRMKVBHD",base))
}

ambig_minus <- function(ambig,ref){ # http://www.virology.wisc.edu/acp/CommonRes/SingleLetterCode.html
    ambig <- toupper(ambig)
    ref   <- toupper(ref)
    if(ambig=="S"){
        if(ref=="G") return("C")
        else if(ref=="C") return("G")
        else return(ambig)
    }else if(ambig=="W"){
        if(ref=="A") return("T")
        else if(ref=="T") return("A")
        else return(ambig)
    }else if(ambig=="R"){
        if(ref=="A") return("G")
        else if(ref=="G") return("A")
        else return(ambig)
    }else if(ambig=="Y"){
        if(ref=="C") return("T")
        else if(ref=="T") return("C")
        else return(ambig)
    }else if(ambig=="K"){
        if(ref=="G") return("T")
        else if(ref=="T") return("G")
        else return(ambig)
    }else if(ambig=="M"){
        if(ref=="A") return("C")
        else if(ref=="C") return("A")
        else return(ambig)
    }else return(ambig)
}

retranslate <- function(calls){

    # USER SAMPLE
    coding <- calls[ord_in_cod>0,list(as.numeric(coding_seq),codon,ord_in_cod,user_sample,reference)]
    setnames(coding,"V1","coding_seq")
    push = 0
    #get missing bases for the first frame if incomplete
    while(coding[1,ord_in_cod]!=1){
        coding<-rbind(coding,cod_table[coding_seq==(as.numeric(coding[1,coding_seq])-1),list(coding_seq=as.numeric(coding_seq),codon,ord_in_cod,user_sample=seq,reference=seq)])
        setkey(coding,coding_seq)
        push = push +1
    }
    #get missing bases for the last frame if incomplete
    while(coding[nrow(coding),ord_in_cod] != 3){
        coding<-rbind(coding,cod_table[coding_seq==(as.numeric(coding[nrow(coding),coding_seq])+1),list(coding_seq=as.numeric(coding_seq),codon,ord_in_cod,user_sample=seq,reference=seq)])
        setkey(coding,coding_seq)
    }
    coding[,user_sample:=ambig_minus(ambig=user_sample,ref=reference),by=1:nrow(coding)]
    trans <- translate(coding[user_sample != '-',user_sample],frame = (coding[1,ord_in_cod]-1), NAstring = "X", ambiguous = F)
    #Shift annotation of codons by '-'s
    ord_sample<-rep(c(1,2,3),length(trans))
    suppressWarnings(trans_sample <- aaa(trans))
#! # calls[ord_in_cod>0,aa_sample := rep(trans,each=3)[(1+push):(length(aa_sample)+push)]]

    # USER MUT
    coding <- calls[ord_in_cod>0,list(as.numeric(coding_seq),codon,ord_in_cod,user_mut,reference)]
    setnames(coding,"V1","coding_seq")
    push = 0
    while(coding[1,ord_in_cod]!=1){
        coding<-rbind(coding,cod_table[coding_seq==(as.numeric(coding[1,coding_seq])-1),list(coding_seq=as.numeric(coding_seq),codon,ord_in_cod,user_mut=seq,reference=seq)])
        setkey(coding,coding_seq)
        push = push +1
    }
    while(coding[nrow(coding),ord_in_cod] != 3){
        coding<-rbind(coding,cod_table[coding_seq==(as.numeric(coding[nrow(coding),coding_seq])+1),list(coding_seq=as.numeric(coding_seq),codon,ord_in_cod,user_mut=seq,reference=seq)])
        setkey(coding,coding_seq)
    }
    coding[,user_mut:=ambig_minus(ambig=user_mut,ref=reference),by=1:nrow(coding)]
#! # ord<-rep(c(1,2,3),length(trans))
    trans <- translate(coding[user_mut != '-',user_mut],frame = (coding[1,ord_in_cod]-1), NAstring = "X", ambiguous = F)
    ord_mut<-rep(c(1,2,3),length(trans))
    suppressWarnings(trans_mut <- aaa(trans))

    #create a string as long as trans rep(123),add them to ord in cod where user_mut != "-"s
    calls[ord_in_cod>0        & user_sample == "-", sample_ord_in_cod := 0  ]
    calls[ord_in_cod>0        & user_mut    == "-", mut_ord_in_cod    := 0  ]
    calls[ord_in_cod>0        & user_sample == "-", aa_sample         := "-"]
    calls[ord_in_cod>0        & user_mut    == "-", aa_mut            := "-"]
    calls[ord_in_cod>0        & user_sample != "-", sample_ord_in_cod := ord_sample[(1+push):(length(aa_sample)+push)]]
    calls[ord_in_cod>0        & user_mut    != "-", mut_ord_in_cod    := ord_mut[(1+push):(length(aa_mut)      +push)]]
    calls[sample_ord_in_cod>0 & user_sample != "-", aa_sample         := rep(trans_sample,each=3)[(1+push):(length(aa_sample)+push)]]
    calls[mut_ord_in_cod>0    & user_mut    != "-", aa_mut            := rep(trans_mut,each=3)[(1+push):(length(aa_mut)      +push)]]
#     if(!is.na(g_hetero_del_tab[1])) dels <- as.vector(unlist(apply(g_hetero_del_tab,1,function(x) x[1]:x[2])))
#     else dels <- numeric()
#     if(!is.na(g_hetero_ins_tab[1])) ins <- as.vector(unlist(apply(g_hetero_ins_tab,1,function(x) x[1]:x[2])))
#     else ins <- numeric()
#     if(max(length(dels),length(ins)) != 0){
#         calls[,aa_mut := incorporate_single_vec(calls[["aa_mut"]],ins,dels,"char",T)]
#     }

    return(calls)
}

get_choices <- function(calls){
    choices <- calls[user_sample != "N" & (user_sample != reference | user_mut != reference) & trace_peak != "NA" & !is.na(gen_coord)]
    if (nrow(choices) > 0) {
        #choices <- choices[,`:=` (user_sample=ambig_minus(user_sample,reference),user_mut=ambig_minus(user_mut,reference)),by=1:nrow(choices)]
        choices[,ids:=NA,by=1:nrow(choices)]
        choices[,user_sample:=ambig_minus(user_sample,reference),by=1:nrow(choices)]
        choices[,user_mut:=ambig_minus(user_mut,reference),by=1:nrow(choices)]
        choices[,sample_peak_pct := mround(sample_peak_pct,1),by=1:nrow(choices)]
        choices[,mut_peak_pct := mround(mut_peak_pct,1),by=1:nrow(choices)]

        choices[,`:=`(coding = paste0("c.", coding_seq), protein="")]
    # 1st variant
        # mismatch
        choices[user_sample != reference & user_sample != "-" & reference   != "-", coding := paste0(coding, reference, ">", user_sample)]#, "(", sample_peak_pct, "%)")]
        # ins
        choices[user_sample != reference & reference   == "-",                      coding := paste0(coding, "ins", user_sample)]#,          "(", sample_peak_pct, "%)")]
        # del
        choices[user_sample != reference & user_sample == "-",                      coding := paste0(coding, "del", reference)]#,            "(", sample_peak_pct, "%)")]
        # protein
        choices[aa_sample   != aa_ref,                                              protein:= paste0("p.", aa_ref, codon, aa_sample)]#,      "(", sample_peak_pct, "%)")]
    # 2nd variant
        # mismatch without 1st variant
        choices[user_mut != reference & user_sample == reference   & user_mut    != "-" & reference   != "-",      coding := paste0(coding, reference, ">", user_mut)]#, "(", mut_peak_pct, "%)")]
        # mismatch with 1st variant
        choices[user_mut != reference & user_sample != reference   & user_mut    != user_sample & user_mut != "-", coding := paste0(coding, ">", user_mut)]#,                 "(", mut_peak_pct, "%)")]
        # ins
        choices[user_mut != reference & user_mut    != user_sample & reference   == "-",                           coding := paste0(coding, "ins", user_mut)]#,         "(", mut_peak_pct, "%)")]
        # del
        choices[user_mut != reference & user_mut    != user_sample & user_mut    == "-",                           coding := paste0(coding, "del", reference)]#,         "(", mut_peak_pct, "%)")]
        # protein without 1st variant
        choices[aa_mut   != aa_ref    & aa_sample   == aa_ref,                                                     protein:= paste0("p.", aa_ref, codon, aa_mut)]#,      "(", mut_peak_pct, "%)")]
        # protein with 1st variant
        choices[aa_mut   != aa_ref    & aa_sample   != aa_ref      & aa_sample   != aa_mut,                        protein:= paste0(protein, aa_mut)]#,                  "(", mut_peak_pct, "%)")]
        #choices[aa_mut != aa_ref & aa_sample== aa_ref & aa_mut=="-",                                           protein:= paste0("p.",aa_ref,codon,"fs")]
    }
    setkey(choices,id)
    return(choices)
}

#remove consecutive single base deletions and replace them with one long deletion in table
get_view<-function(calls,choices){

    computeConsecutives <- function(ids){
        ids <- round(ids * 100)
        res <- rle(ids - c(Inf,ids[-length(ids)]))
        ccc <- cumsum(res$lengths)
        ccc[which(res$values != 1 & res$values != 100)] <- 0
        res <- rep(ccc,res$lengths)
        res2 <- res - c(res[-1],0)
        res2[res2 >= 0] <- 0
        return(res - res2)
    }

    squeeze_indels <- function(tab){
        if(nrow(tab) > 0){
            coord <- gsub("c\\.(\\d*).*","\\1",tab$coding)
            nucs <- gsub("c\\.\\d*...(.)","\\1",tab$coding)
            type <- gsub("c\\.\\d*(...).*","\\1",tab$coding)[1]

            if(max(tab$gen_coord) == min(tab$gen_coord)) gen_coord <- paste0(max(tab$gen_coord) + 1,"_",as.numeric(max(tab$gen_coord)))
            else gen_coord <- paste0(max(tab$gen_coord),"_",min(tab$gen_coord))

            if(max(coord) == min(coord)) coding <- paste0("c.",as.numeric(min(coord)) - 1,"_",min(coord),type,ifelse(nrow(tab) > 10,paste0(nrow(tab),"nt"), paste(nucs,collapse = "") ))
            else coding <- paste0("c.",min(coord),"_",max(coord),type, paste(nucs,collapse = ""))

            return(list(id = floor(min(tab$id)),gen_coord = gen_coord,coding = coding,protein = tab$protein[1]))
        } else {
            return(tab)
        }
    }

    choices[,consecutives := computeConsecutives(id) ][,mut_type := gsub("c\\.\\d*(...).*","\\1",coding)]
    indel_tab <- choices[intersect(grep("del|ins",mut_type),which(consecutives != 0)), squeeze_indels(.SD),by = c("mut_type","consecutives")]
    #represent consecutive indels on one line
    if(nrow(indel_tab) > 0){
        choices <- choices[union(grep("del|ins",mut_type,invert = T),which(consecutives == 0)),]
        choices <- rbind(choices,indel_tab,fill=TRUE)
    }
    #identify frame shifts and inframe indels
    for(i in grep("del|ins",choices$mut_type)){
        seq <- gsub("c\\.\\d*_*\\d*...(.)","\\1",choices[i,]$coding)
        #if(length(seq) > 10,paste0(length(seq),"nt")
        if((str_length(seq) %% 3)!=0){
            prot <- gsub("(p\\....\\d*).*","\\1",choices[i]$protein)
            aa <- gsub("p\\.(...)\\d*.*","\\1",choices[i]$protein)
            cod <- as.numeric(gsub("p\\....(\\d*).*","\\1",choices[i]$protein))
            while((calls[codon == cod][1]$aa_ref == calls[codon == cod][1]$aa_mut)&
                  (calls[codon == cod][1]$aa_ref == calls[codon == cod][1]$aa_sample)&
                  (!is.na(calls[codon == cod][1]$aa_sample))) {cod = cod +1}
            choices[i,]$protein = paste0("p.",aa,cod, "fs")

        }else{ #in frame
            if(choices[i,]$mut_type == "ins"){
                from <- as.numeric(calls[choices[i]$id]$codon)
                to   <- as.numeric(calls[choices[i]$id]$codon) +1
                choices[i,]$protein = paste0("p.",calls[codon==from,][1]$aa_ref,from,"_",
                                             calls[codon==to,][1]$aa_ref,to,choices[i,]$mut_type,
                                             paste(translate(strsplit(seq,"")[[1]]),collapse = ""))
            }
            if(choices[i,]$mut_type=="del"){

                from <- as.numeric(calls[choices[i]$id]$codon)
                to   <- as.numeric(calls[choices[i]$id]$codon) + nchar(seq)/3 -1

                #from <- as.numeric(g_calls[choices[i]$id]$codon) - 10
                #to   <- as.numeric(g_calls[choices[i]$id]$codon) + nchar(seq)/3 + 10
                #lapply(g_calls[codon %in% c(from:to) & ord_in_cod == 1]$aa_ref,mya)

                choices[i,]$protein = paste0("p.",calls[codon==from,][1]$aa_ref,from,"_",
                                             calls[codon==to,][1]$aa_ref,to,choices[i,]$mut_type)
            }
        }
    }
    #identify duplications (special kind of insertions)
    for(i in grep("ins",choices$mut_type)){
        seq <- gsub("c\\.\\d*_*\\d*...(.)","\\1",choices[i,]$coding)
        if(floor(choices[i,]$id) == choices[i,]$id) prev_seq <- paste0(calls[-(nchar(seq) - 1):0 + choices[i,]$id - 1,]$reference,collapse = "")
        else prev_seq <- paste0(calls[-(nchar(seq) - 1):0 + choices[i,]$id,]$reference,collapse = "")
        if(seq == prev_seq) {
            choices[i,coding := gsub("ins","dup",coding)]
            #the coordinates are changed to the sequence that is duplicated #! find teste for these
            if(nchar(seq)>1){
                choices[i,]$coding <- paste0("c.",calls[floor(choices[i,]$id)-nchar(seq)+1]$coding_seq,"_",calls[floor(choices[i,]$id),]$coding_seq,"dup",seq)
            }
        }else{
            if(nchar(seq)==1)
                choices[i,]$coding <- paste0("c.",calls[floor(choices[i]$id),]$coding_seq,"_",calls[ceiling(choices[i]$id),]$coding_seq,"ins",seq)
        }
    }

    setkey(choices,id)
    return(choices)
}
mya <- function(x){
    if(is.na(x)){return(".")}
    else{ if(x == "-"){return("-")}
        else{return(a(x))}
    }
}
mround <- function(x,base){
    base*round(x/base)
}

report_hetero_indels <- function(calls){
    rev <- !is.null(calls[["call_rev"]])
    if(rev) {
        secondary_seq <- gsub("[ -]","",paste(get_consensus_mut(calls[["mut_call_fwd"]],calls[["mut_call_rev"]],calls[,list(iA_fwd,iC_fwd,iG_fwd,iT_fwd,iA_rev,iC_rev,iG_rev,iT_rev)],calls[["user_sample"]]),collapse = ""))
    } else secondary_seq <- gsub("[ -]","",paste(calls[["mut_call_fwd"]],collapse = ""))
    primary_seq <- gsub("[ -]","",paste(calls[["user_sample"]],collapse = ""))
    hetero_indel_aln <- pairwiseAlignment(primary_seq, secondary_seq,type = "overlap",substitutionMatrix = sm,gapOpening = -15, gapExtension = -1)

    hetero_ins_tab <- stringi::stri_locate_all_regex(compareStrings(hetero_indel_aln),"[\\-]+")[[1]] + start(pattern(hetero_indel_aln)) - 1
    hetero_del_tab <- stringi::stri_locate_all_regex(compareStrings(hetero_indel_aln),"[\\+]+")[[1]] + start(pattern(hetero_indel_aln)) - 1

    is.in.primery <- apply(hetero_ins_tab,1,function(x) all(x %in% which(calls[["user_sample"]] == "-")))

    if(length(hetero_ins_tab[which(!is.in.primery),2]) > 0){
        move_vec <- numeric(nchar(primary_seq))
        move_vec[hetero_ins_tab[which(!is.in.primery),2] + 1] <- -(hetero_ins_tab[which(!is.in.primery),2]-hetero_ins_tab[which(!is.in.primery),1])-1
        move_vec <- cumsum(move_vec)
        hetero_del_tab <- apply(hetero_del_tab,c(1,2),function(x) x + move_vec[x])
        hetero_ins_tab <- apply(hetero_ins_tab,c(1,2),function(x) x + move_vec[x])
    }

    is.in.reference <- apply(hetero_del_tab,1,function(x) all(x %in% which(calls[["reference"]] == "-")))

    ins_counts <- sum(hetero_ins_tab[which(!is.in.primery),2]-hetero_ins_tab[which(!is.in.primery),1]+1,na.rm = T) + sum(hetero_del_tab[which(is.in.reference),2]-hetero_del_tab[which(is.in.reference),1]+1,na.rm = T)
    del_counts <- sum(hetero_ins_tab[which(is.in.primery),2]-hetero_ins_tab[which(is.in.primery),1]+1,na.rm = T) + sum(hetero_del_tab[which(!is.in.reference),2]-hetero_del_tab[which(!is.in.reference),1]+1,na.rm = T)
    # if(nrow(hetero_ins_tab) > 0) g_minor_het_insertions <<- data.table::data.table(pos = )
    if(nrow(hetero_ins_tab) > 0) {
        offset <- -start(pattern(hetero_indel_aln)) + start(subject(hetero_indel_aln))
        minor_het_insertions <<- data.table::data.table(pos = hetero_ins_tab[which(!is.in.primery),1],seq = stri_sub(secondary_seq,hetero_ins_tab[which(!is.in.primery),1] + offset,hetero_ins_tab[which(!is.in.primery),2] + offset))
    }
        else minor_het_insertions <<- data.table::data.table()
    hetero_indel_aln <<- hetero_indel_aln
    hetero_indel_pid <<- round(pid(hetero_indel_aln),1)
    hetero_ins_tab   <<- hetero_ins_tab
    hetero_del_tab   <<- hetero_del_tab
    hetero_indel_report <<- paste0("alignment % id : ",hetero_indel_pid,"%\nins/del counts : ",ins_counts," / ",del_counts)
    if((ins_counts > 0)||(del_counts > 0)){
        indels_present <- TRUE
    }else{
        indels_present <- FALSE
    }
    return(list(indels_present=indels_present,minor_het_insertions=minor_het_insertions,hetero_indel_aln=hetero_indel_aln,hetero_ins_tab=hetero_ins_tab,hetero_del_tab=hetero_del_tab,hetero_indel_pid=hetero_indel_pid,hetero_indel_report=hetero_indel_report))
}

get_consensus_mut <- function(mut_fwd,mut_rev,intens_tab,primery_seq){
#     if(length(which(primery_seq == "-")) > 0){
#         mut_fwd <- mut_fwd[-which(primery_seq == "-")]
#         mut_rev <- mut_rev[-which(primery_seq == "-")]
#     }
    names <- structure(c(1,1:4),names = c("N","A","C","G","T"))
    pa <- pairwiseAlignment(gsub("[ -]","",paste(mut_fwd,collapse = "")), gsub("[ -]","",paste(mut_rev,collapse = "")),type = "overlap",substitutionMatrix = sm,gapOpening = -10, gapExtension = -1)
    fwd <- strsplit(as.character(pattern(pa)),"")[[1]]
    rev <- strsplit(as.character(subject(pa)),"")[[1]]
    fwd_i <- numeric(length(fwd))
    rev_i <- numeric(length(rev))
    fwd_start <- min(which(mut_fwd != "-")) + start(pattern(pa)) - 1
    rev_start <- min(which(mut_rev != "-")) + start(subject(pa)) - 1
    fwd_matrix <- as.matrix(intens_tab[,1:4,with = F])
    rev_matrix <- as.matrix(intens_tab[,5:8,with = F])
    fwd_i[which(fwd != "-")] <- diag(fwd_matrix[start(pattern(pa)):nrow(fwd_matrix),names[fwd[which(fwd != "-")]]])
    rev_i[which(rev != "-")] <- diag(rev_matrix[start(subject(pa)):nrow(rev_matrix),names[rev[which(rev != "-")]]])
    cons <- fwd
    cons[which(rev_i > fwd_i)] <- rev[which(rev_i > fwd_i)]
    return(cons)
}

#reconstructs user_mut and mut_peak_pct by shifting user_call_fwd and user_call_rev by detected indels indels
incorporate_hetero_indels_func <- function(calls,hetero_del_tab,hetero_ins_tab,minor_het_insertions){
    if(!is.na(hetero_del_tab[1])) dels <- as.vector(unlist(apply(hetero_del_tab,1,function(x) x[1]:x[2])))
    else dels <- numeric()
    if(!is.na(hetero_ins_tab[1])) ins <- as.vector(unlist(apply(hetero_ins_tab,1,function(x) x[1]:x[2])))
    else ins <- numeric()
    if(max(length(dels),length(ins)) != 0){
        calls[,het_mut_call_fwd     := incorporate_single_vec(calls[["mut_call_fwd"]],ins,dels,"char",T,calls[["sample_peak_base_fwd"]])]
        calls[,het_mut_peak_pct_fwd := incorporate_single_vec(calls[["mut_peak_pct_fwd"]],ins,dels,"num",T)]
#         calls[,het_mut_s2n_abs_fwd := incorporate_single_vec(calls[["mut_s2n_abs_fwd"]],ins,dels,"num",T)]
        calls[set_by_user == FALSE, c("user_mut","mut_peak_pct") := list(het_mut_call_fwd,het_mut_peak_pct_fwd)]
        if(any(colnames(calls) == "call_rev")){
            calls[,het_mut_call_rev     := incorporate_single_vec(calls[["mut_call_rev"]],ins,dels,"char",F,calls[["sample_peak_base_rev"]])]
            calls[,het_mut_peak_pct_rev := incorporate_single_vec(calls[["mut_peak_pct_rev"]],ins,dels,"num",F)]
#             calls[,het_mut_s2n_abs_rev := incorporate_single_vec(calls[["mut_s2n_abs_rev"]],ins,dels,"num",F)]
            calls[(quality_fwd < quality_rev | user_mut == "-") & set_by_user == FALSE, c("user_mut","mut_peak_pct") := list(het_mut_call_rev,het_mut_peak_pct_rev)]
        }
    }
    if(nrow(minor_het_insertions[!is.na(pos)]) > 0){
        get_ins_data_table <- function(pos,seq){
            ins_seq <- strsplit(seq,"")[[1]]
            ins_tab <- calls[rep(pos-1,length(ins_seq)),]
            ins_tab[,id := id + seq_along(id)/100]
            ins_tab[,user_sample := "-"][,reference := "-"][,user_mut := ins_seq]
            ins_tab[,`:=`(iA_fwd=0,iC_fwd=0,iG_fwd=0,iT_fwd=0,ord_in_cod=4)]
            if("call_rev" %in% row.names(calls)){
                ins_tab[,`:=`(iA_rev=0,iC_rev=0,iG_rev=0,iT_rev=0)]
            }
            return(ins_tab)
        }

        if(nrow(minor_het_insertions[!is.na(pos)])>0){
            ins_tabs <- lapply(1:nrow(minor_het_insertions[!is.na(pos),]),function(x) get_ins_data_table(minor_het_insertions[!is.na(pos),]$pos[x],minor_het_insertions$seq[x]))
            minor_het_insertions$added <<- lapply(1:nrow(minor_het_insertions[!is.na(pos),]),function(x) paste0(ins_tabs[[x]]$id,collapse= " "))
            #g_minor_het_insertions$added = rbindlist(ins_tabs)$id;
            calls <- rbindlist(c(list(calls),ins_tabs))
        }
    }
    return(calls=calls,minor_het_insertions=minor_het_insertions)
}

incorporate_single_vec <- function(vec,ins,dels,type,fwd,primarySeq){
    orig_vec <- vec
    if(type == "num") elem <- 0
    else elem <- "-"
    new_vec <- rep(elem,length(vec))
    if(length(ins) > 0) {
        if(fwd) {
            vec <- vec[-ins]
            vec <- c(vec,rep(elem,length(ins)))
        } else {
            vec <- vec[-(ins - length(ins))]
            vec <- c(rep(elem,length(ins)),vec)
        }
    }
    if(fwd) {
        vec <- vec[1:min(length(vec),length(new_vec) - length(dels))]
        new_vec[setdiff(seq_along(new_vec),dels)] <- vec
    } else {
        vec <- vec[1 + length(vec) -  min(length(vec),length(new_vec) - length(dels)):1]
        new_vec[setdiff(seq_along(new_vec),dels)] <- vec
    }
    if(type == "char" && !is.null(primarySeq)){
        move_vec <- numeric(length(new_vec))
        if(fwd){
            move_vec[ins] <- -1
            move_vec[dels] <- 1
            move_vec <- cumsum(move_vec)
            move_vec[which(rep(rle(move_vec)$length,rle(move_vec)$length) == 1)] <- move_vec[which(rep(rle(move_vec)$length,rle(move_vec)$length) == 1) + 1]
        } else {
            move_vec[ins] <- 1
            move_vec[dels] <- -1
            move_vec <- rev(cumsum(rev(move_vec)))
            move_vec[which(rep(rle(move_vec)$length,rle(move_vec)$length) == 1)] <- move_vec[which(rep(rle(move_vec)$length,rle(move_vec)$length) == 1) - 1]
        }
        replace <- setdiff(which(orig_vec == primarySeq) + move_vec[which(orig_vec == primarySeq)],c(ins,dels))
        replace <- replace[replace < length(primarySeq)]
        replace <- replace[replace > 0]
        new_vec[replace] <- primarySeq[replace]
    }
    return(new_vec)
}

add_intensities <- function(added){
    #update intensities
    id <- g_calls[id == as.integer(added[1]),]$trace_peak + 6
    add<-data.table("id"=id + (1:(length(added)*12)/1000),"A"=0,"C"=0,"G"=0,"T"=0)
    g_intens     <<- rbind(g_intens,     add)
    setkey(g_intens,     id)
    if("call_rev" %in% colnames(g_calls)){
        g_intens_rev <<- rbind(g_intens_rev, add)
        setkey(g_intens_rev, id)

    }
    #intens_rev must match intens (hopefully they do otherwise its a bigger problem)
    #update peak positions in calls table
    g_calls$trace_peak<<-seq(from = g_calls[1]$trace_peak, by = 12, length.out = nrow(g_calls))
    g_calls$trace_peak_rev<<-seq(from = g_calls[1]$trace_peak, by = 12, length.out = nrow(g_calls))
    #update intrex
    g_intrexdat$intrex     <- setnames(g_calls[!is.na(exon_intron),list(max(id)-min(id)+1,min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","length","trace_peak","end"))
    g_intrexdat$intrex     <- setnames(merge(g_intrexdat$intrex,g_calls[,list(id,trace_peak)],by="trace_peak"),"trace_peak","start")
    g_intrexdat       <<- splice_variants(g_intrexdat)
    g_intrexdat$max_x <<- nrow(g_intens)
    return(paste0(add$id,collapse= " "))
}

remove_intensities <- function(added){
    #update intensities (this operation takes too long)
    g_intens <<- g_intens[!id %in% as.numeric(str_split(g_minor_het_insertions$ins_added," ")[[1]]),]
    #update peak positions in calls table
    g_calls$trace_peak<<-seq(from = g_calls[1]$trace_peak, by = 12, length.out = nrow(g_calls))
    if("call_rev" %in% colnames(g_calls)){
        g_intens_rev <<-g_intens_rev[! id %in% as.numeric(str_split(g_minor_het_insertions$ins_added," ")[[1]]),]
        g_calls$trace_peak_rev<<-seq(from = g_calls[1]$trace_peak, by = 12, length.out = nrow(g_calls))
    }
    #update intrex
    g_intrexdat$intrex     <- setnames(g_calls[!is.na(exon_intron),list(max(id)-min(id)+1,min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","length","trace_peak","end"))
    g_intrexdat$intrex     <- setnames(merge(g_intrexdat$intrex,g_calls[,list(id,trace_peak)],by="trace_peak"),"trace_peak","start")
    g_intrexdat       <<- splice_variants(g_intrexdat)
    g_intrexdat$max_x <<- nrow(g_intens)
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
         if (mut_peak == 0 ) return(list("-",mut_peak))
    else if (mut_peak == iA) return(list("A",mut_peak))
    else if (mut_peak == iC) return(list("C",mut_peak))
    else if (mut_peak == iG) return(list("G",mut_peak))
    else if (mut_peak == iT) return(list("T",mut_peak))
    else                     return(list("-",mut_peak))
}

# i_w_p <- function(p,iA,iC,iG,iT,pA,pC,pG,pT){
#     mut_peak <- (sort(c(iA,iC,iG,iT),decreasing = TRUE)[p])
#          if (mut_peak == 0 ) return(list(" ",mut_peak,pA))
#     else if (mut_peak == iA) return(list("A",mut_peak,pA))
#     else if (mut_peak == iC) return(list("C",mut_peak,pC))
#     else if (mut_peak == iG) return(list("G",mut_peak,pG))
#     else if (mut_peak == iT) return(list("T",mut_peak,pT))
# }


get_expected_het_indels <- function(calls){
    min_het_pct <- 0.04

    rev <- !is.null(calls[["call_rev"]])

    if(rev) intens_tab <- calls[,list(iA_fwd,iC_fwd,iG_fwd,iT_fwd,iA_rev,iC_rev,iG_rev,iT_rev)]
    else intens_tab <- calls[,list(iA_fwd,iC_fwd,iG_fwd,iT_fwd)]

    fwd_matrix <- as.matrix(intens_tab[,1:4,with = F])
    fwd_matrix2 <- fwd_matrix / rowSums(fwd_matrix)
    fwd_matrix2[which(is.na(fwd_matrix2))] <- 0
    sec_fwd_vec <- apply(fwd_matrix2,1,function(x) max(x[-which.max(x)]))

    if(rev){
        rev_matrix <- as.matrix(intens_tab[,5:8,with = F])
        rev_matrix2 <- rev_matrix / rowSums(rev_matrix)
        rev_matrix2[which(is.na(rev_matrix2))] <- 0
        sec_rev_vec <- apply(rev_matrix2,1,function(x) max(x[-which.max(x)]))
        sec_vec <- c(sec_fwd_vec,sec_rev_vec)
    } else {
        sec_vec <- sec_fwd_vec
    }

    hst <- hist(sec_vec,breaks = 100,plot = F)

    xz <- as.zoo(hst$density)
    min <- rollapply(xz, 9, function(x) which.min(x)==5)
    max <- rollapply(xz, 9, function(x) which.max(x)==5)

    max <- index(max)[which(max)]
    max <- max[max > which(hst$breaks > min_het_pct)[1]]
    max_dens <- rollapply(xz, 9, sum)[max]
    best_max <- index(max_dens)[which.max(max_dens)]

    min <- index(min)[which(min)]
    min <- min[min < best_max]
    if(length(min) == 0){
        g_expected_het_indel <- list(min = 0,max = 0,hist = hst)
        return(g_expected_het_indel)
    }

    best_min <- rev(min)[which.min(rev(hst$density[min]))]

    if(max_dens[which(index(max_dens) == best_max)] > 8){
        g_expected_het_indel <- list(min = hst$breaks[best_min + 1],max = hst$breaks[best_max + 1],hist = hst)
        return(g_expected_het_indel)
    } else  {
        g_expected_het_indel <- list(min = 0,max = 0,hist = hst)
        return(g_expected_het_indel)
    }
    
}
