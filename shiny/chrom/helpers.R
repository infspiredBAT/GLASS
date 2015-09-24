annotate_calls <- function(calls,intens,intens_rev){

    #contains codons table
    load("../../data/codons.rdata")
    calls <- merge(x = calls, y = cod_table[,list(gen_coord,codon,ord_in_cod,coding_seq,AA_ref=AA,AA_mod=AA)], by = "gen_coord", all.x = TRUE)
    cod_table <<- cod_table
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
    # calls[,mut_call_fwd:=sample_peak_base_fwd]
    calls[,mut_call_fwd:=call]
    if("call_rev" %in% colnames(calls)){
        calls[,c("sample_peak_base_rev","sample_peak_abs_rev") := i_wo_p(1,iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]
        calls[,sample_peak_pct_rev := ((100/ref_peak_abs_rev)*sample_peak_abs_rev)]
        calls[,c("mut_peak_base_rev","mut_peak_abs_rev") := i_wo_p(2,iA_rev,iC_rev,iG_rev,iT_rev),by=1:nrow(calls)]
        calls[,mut_peak_pct_rev := ((100/ref_peak_abs_rev)*mut_peak_abs_rev)]
        calls[,mut_s2n_abs_rev:=mut_peak_abs_rev/noise_abs_rev]
        # calls[,mut_call_rev:=sample_peak_base_rev]
        calls[,mut_call_rev:=call_rev]

         calls[,sample_peak_pct := max(c(sample_peak_pct_fwd,sample_peak_pct_rev),na.rm=TRUE),by=1:nrow(calls)]
#         calls[,mut_peak_pct := mean(c(mut_peak_pct_fwd,mut_peak_pct_rev),na.rm=TRUE),by=1:nrow(calls)]
    } else {
         calls[,sample_peak_pct := sample_peak_pct_fwd]
#         calls[,mut_peak_pct := mut_peak_pct_fwd]
    }
    calls[,set_by_user:=FALSE]
    calls[,cons:=user_sample]
    return(calls)
}

call_variants <- function(calls, qual_thres, mut_min, s2n_min){
    # reset all but set_by_user
    # calls[set_by_user == FALSE, mut_call_fwd := sample_peak_base_fwd]
    calls[set_by_user == FALSE, mut_call_fwd := call]
    calls[set_by_user == FALSE, user_sample := cons]
    calls[set_by_user == FALSE, user_mut := cons]

    calls[rm7qual <= qual_thres & set_by_user == FALSE, user_sample := "N"]

    # mut
    if("call_rev" %in% colnames(calls)) {
        # reset all but set_by_user
        # calls[set_by_user == FALSE, mut_call_rev := sample_peak_base_rev]
        calls[set_by_user == FALSE, mut_call_rev := call_rev]
        calls[
            mut_peak_pct_fwd >= mut_min
            & mut_s2n_abs_fwd >= s2n_min
            ,mut_call_fwd := mut_peak_base_fwd
            ]
        calls[
            mut_peak_pct_rev >= mut_min
            & mut_s2n_abs_rev >= s2n_min
            ,mut_call_rev := mut_peak_base_rev
            ]
        calls[
            set_by_user == FALSE
            & mut_call_fwd != call
            ,c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
            ]
        calls[
            ((
            mut_peak_pct_rev > mut_peak_pct_fwd
            & mut_s2n_abs_rev > mut_s2n_abs_fwd
            )
            | user_mut == "-"
            )
            & set_by_user == FALSE
            & mut_call_rev != call_rev
            ,c("user_mut","mut_peak_pct") := list(mut_call_rev,mut_peak_pct_rev)
            ]
    } else {
        calls[
            mut_peak_pct_fwd >= mut_min
            & mut_s2n_abs_fwd >= s2n_min
            ,mut_call_fwd := mut_peak_base_fwd
            ]
        calls[
            set_by_user == FALSE
            ,user_mut := mut_call_fwd
            ]
    }
    return(calls)
}

complement <- function(base){
    return (chartr("ATGC","TACG",base))
}
retranslate <- function(calls){
  coding <- calls[ord_in_cod>0,list(coding_seq,codon,ord_in_cod,user_sample)]
  push = 0
  while(coding[1,ord_in_cod]!=1){
    coding<-rbind(coding,cod_table[coding_seq==(as.numeric(coding[1,coding_seq])-1),list(coding_seq,codon,ord_in_cod,user_sample=seq)])
    setkey(coding,coding_seq)
    push = push +1
  }
  while(coding[nrow(coding),ord_in_cod] != 3){
      coding<-rbind(coding,cod_table[coding_seq==(as.numeric(coding[nrow(coding),coding_seq])+1),list(coding_seq,codon,ord_in_cod,user_sample=seq)])
      setkey(coding,coding_seq)
  }
  trans <- translate(coding[,user_sample],frame = (coding[1,ord_in_cod]-1))
  return(calls[ord_in_cod>0,AA_mod := rep(trans,each=3)[(1+push):(length(AA_mod)+push)]])
}

get_choices <- function(calls){
    choices <- calls[user_sample != "N" & (user_sample != reference | user_mut != reference) & trace_peak != "NA" & !is.na(gen_coord)]

    choices[user_sample == "-",user_sample := "del"]
    choices[user_mut == "-",user_mut := "del"]

    choices[,sample_peak_pct:=mround(sample_peak_pct,5),by=1:nrow(choices)]
    choices[,mut_peak_pct:=mround(mut_peak_pct,5),by=1:nrow(choices)]

    choices[,`:=`(coding=paste0("c.",coding_seq,reference,">"),protein="")]
    choices[user_sample != reference,coding:=paste0(coding,user_sample,"(",sample_peak_pct,"%)")]
    choices[user_mut != reference & user_mut != user_sample,coding:=paste0(coding,user_mut,"(",mut_peak_pct,"%)")]
    choices[AA_ref!="",protein:=paste0("p.",codon,toupper(AA_ref),">",AA_mod)]

    return(choices)
}

mround <- function(x,base){
    base*round(x/base)
}

report_hetero_indels <- function(calls){
    rev <- !is.null(calls[["call_rev"]])
    if(rev) {
        secondary_seq <- paste(get_consensus_mut(calls[["mut_call_fwd"]],calls[["mut_call_rev"]],calls[,list(iA_fwd,iC_fwd,iG_fwd,iT_fwd,iA_rev,iC_rev,iG_rev,iT_rev)]),collapse = "")
    } else secondary_seq <- gsub("[ -]","",paste(calls[["mut_call_fwd"]],collapse = ""))
    primary_seq <- gsub("[ -]","",paste(calls[["user_sample"]],collapse = ""))
    hetero_indel_aln <- pairwiseAlignment(primary_seq, secondary_seq,type = "local",substitutionMatrix = sm,gapOpening = -20, gapExtension = -0.5)
    # writePairwiseAlignments(indel_align, block.width = 150)

    # cat("identified deletions:\n")
    hetero_del_tab <- stringi::stri_locate_all_regex(compareStrings(hetero_indel_aln),"[\\+]+")[[1]] + start(pattern(hetero_indel_aln)) - 1
#     if(is.na(hetero_indel_tab[1])) cat("no deletions")
#     else print(hetero_indel_tab)
#     cat("\nidentified insertions:\n")
    hetero_ins_tab <- stringi::stri_locate_all_regex(compareStrings(hetero_indel_aln),"[\\-]+")[[1]] + start(pattern(hetero_indel_aln)) - 1
#     if(is.na(hetero_indel_tab[1])) cat("no insertions")
#     else print(hetero_indel_tab)

    hetero_indel_pid <- round(pid(hetero_indel_aln),1)
    return(list(hetero_indel_aln, hetero_indel_pid, hetero_ins_tab, hetero_del_tab))
}

get_consensus_mut <- function(mut_fwd,mut_rev,intens_tab){
    names <- structure(1:4,names = c("A","C","G","T"))
    pa <- pairwiseAlignment(gsub("[ -]","",paste(mut_fwd,collapse = "")), gsub("[ -]","",paste(mut_rev,collapse = "")),type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
    fwd <- strsplit(as.character(pattern(pa)),"")[[1]]
    rev <- strsplit(as.character(subject(pa)),"")[[1]]
    fwd_i <- numeric(length(fwd))
    rev_i <- numeric(length(rev))
    fwd_start <- min(which(mut_fwd != " ")) + start(pattern(pa)) - 1
    rev_start <- min(which(mut_rev != " ")) + start(subject(pa)) - 1
    fwd_i[which(fwd != "-")] <- diag(as.matrix(intens_tab[,1:4,with = F])[fwd_start:(fwd_start + length(which(fwd != "-"))),names[fwd[which(fwd != "-")]]])
    rev_i[which(rev != "-")] <- diag(as.matrix(intens_tab[,5:8,with = F])[rev_start:(rev_start + length(which(rev != "-"))),names[rev[which(rev != "-")]]])
    cons <- fwd
    cons[which(rev_i > fwd_i)] <- rev[which(rev_i > fwd_i)]
    return(cons)
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
