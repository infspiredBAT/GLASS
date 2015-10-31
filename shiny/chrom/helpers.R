annotate_calls <- function(calls,intens,intens_rev){

    #contains codons table
    load("data/codons.rdata")
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

        suppressWarnings(calls[,sample_peak_pct := max(c(sample_peak_pct_fwd,sample_peak_pct_rev),na.rm=TRUE),by=1:nrow(calls)])
        calls[,mut_peak_pct := mean(c(mut_peak_pct_fwd,mut_peak_pct_rev),na.rm=TRUE),by=1:nrow(calls)]
    } else {
        calls[,sample_peak_pct := sample_peak_pct_fwd]
        calls[,mut_peak_pct := mut_peak_pct_fwd]
    }
    calls[,set_by_user:=FALSE]
    calls[,user_sample_orig:=user_sample]
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

call_variants <- function(calls, qual_thres, mut_min, s2n_min){
    # reset all but set_by_user
    # calls[set_by_user == FALSE, mut_call_fwd := sample_peak_base_fwd]
    calls[set_by_user == FALSE, mut_call_fwd := ambig_minus(call,reference),by=1:nrow(calls[set_by_user==FALSE,])]
    calls[set_by_user == FALSE, user_sample := user_sample_orig]
    calls[set_by_user == FALSE, user_mut := user_sample_orig]

    # mut
    if("call_rev" %in% colnames(calls)) {
        # reset all but set_by_user
        # calls[set_by_user == FALSE, mut_call_rev := sample_peak_base_rev]
        calls[set_by_user == FALSE, mut_call_rev := ambig_minus(call_rev,reference),by=1:nrow(calls[set_by_user==FALSE,])]
        calls[
            mut_peak_pct_fwd >= mut_min
            & mut_s2n_abs_fwd >= s2n_min
            & mut_peak_base_fwd != reference
            #& quality_fwd >= qual_thres
            ,mut_call_fwd := mut_peak_base_fwd
            ]
        calls[
            mut_peak_pct_rev >= mut_min
            & mut_s2n_abs_rev >= s2n_min
            & mut_peak_base_rev != reference
            #& quality_rev >= qual_thres
            ,mut_call_rev := mut_peak_base_rev
            ]
#         calls[
#             set_by_user == FALSE
#             # & mut_call_fwd != call
#             ,c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
#             ]
        calls[mut_call_fwd!=reference
              &mut_call_rev==reference
              &quality_fwd > qual_thres
              &set_by_user == FALSE,
              c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
              ]
        calls[mut_call_fwd==reference
              &mut_call_rev!=reference
              &quality_rev > qual_thres
              &set_by_user == FALSE,
              c("user_mut","mut_peak_pct") := list(mut_call_rev,mut_peak_pct_rev)
              ]

        calls[
            mut_call_fwd!=reference
            &mut_call_rev!=reference
            & set_by_user == FALSE
            # & mut_call_rev != call_rev
            ,c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
            ]
        calls[
            mut_call_fwd!=reference
            &mut_call_rev!=reference
            &quality_rev > quality_fwd
            & set_by_user == FALSE
            # & mut_call_rev != call_rev
            ,c("user_mut","mut_peak_pct") := list(mut_call_rev,mut_peak_pct_rev)
            ]

    } else {
        calls[
            mut_peak_pct_fwd >= mut_min
            & mut_s2n_abs_fwd >= s2n_min
            & mut_peak_base_fwd != reference
            ,mut_call_fwd := mut_peak_base_fwd
            ]
        calls[
            set_by_user == FALSE
            & mut_call_fwd != call
            ,c("user_mut","mut_peak_pct") := list(mut_call_fwd,mut_peak_pct_fwd)
            ]
    }
    calls[quality < qual_thres & set_by_user == FALSE, c("user_sample","user_mut") := "N"]
    return(calls)
}

complement <- function(base){
    return (chartr("ATGC","TACG",base))
}
ambig_minus <- function(ambig,ref){ # http://www.virology.wisc.edu/acp/CommonRes/SingleLetterCode.html
    ambig <- toupper(ambig)
    ref   <- toupper(ref)
    if(ambig=="S"){
        if(ref=="G") return("C")
        if(ref=="C") return("G")
    }else if(ambig=="W"){
        if(ref=="A") return("T")
        if(ref=="T") return("A")
    }else if(ambig=="R"){
        if(ref=="A") return("G")
        if(ref=="G") return("A")
    }else if(ambig=="Y"){
        if(ref=="C") return("T")
        if(ref=="T") return("C")
    }else if(ambig=="K"){
        if(ref=="G") return("T")
        if(ref=="T") return("G")
    }else if(ambig=="M"){
        if(ref=="A") return("C")
        if(ref=="C") return("A")
    }else return(ambig)
}
retranslate <- function(calls){

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

    suppressWarnings(trans <- aaa(trans))
    calls[ord_in_cod>0,aa_sample := rep(trans,each=3)[(1+push):(length(aa_sample)+push)]]

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
    ord<-rep(c(1,2,3),length(trans))
    trans <- translate(coding[user_mut != '-',user_mut],frame = (coding[1,ord_in_cod]-1), NAstring = "X", ambiguous = F)
    suppressWarnings(trans <- aaa(trans))
    #create a string as long as trans rep(123),add them to ord in cod where user_mut != "-"s
    calls[ord_in_cod>0 & user_mut == "-",mut_ord_in_cod:=0]
    calls[ord_in_cod>0 & user_mut == "-",aa_mut:="-"]
    calls[ord_in_cod>0 &user_mut!="-",mut_ord_in_cod := ord[(1+push):(length(aa_mut)+push)]]
    calls[mut_ord_in_cod>0 & user_mut != "-",aa_mut := rep(trans,each=3)[(1+push):(length(aa_mut)+push)]]
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
        choices[,sample_peak_pct := mround(sample_peak_pct,5),by=1:nrow(choices)]
        choices[,mut_peak_pct := mround(mut_peak_pct,5),by=1:nrow(choices)]

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
get_view<-function(choices){
    reps<-rle(choices$user_mut)
    reps$csum <- cumsum(reps$length)
    dtr <- data.table(length=reps$length,values=reps$values,csum=reps$csum)
    if(nrow(dtr[values=="-"&length>=2,])>=1){
        dtr[values=="-"&length>=2,`:=`(start=csum-length+1,end=csum)]
        deletions <- dtr[values=="-"&length>=2,
                         list(id=choices[start,]$id,
                              ids=gsub(",","",toString(choices[,id[start:end]])), #Might need this for locks and resets
                              start=start,                                        #to now which rows I'm replacing
                              end=end,
                              gen_coord=paste0(choices[start,]$gen_coord,"_",
                                               choices[end,]$gen_coord),
                              coding=paste0("c.",choices[start]$coding_seq,"_",
                                            choices[end,]$coding_seq,"del",
                                            gsub(", ","",toString(choices[,reference[start:end]]))),
                              protein=paste0("p.",choices[start,]$aa_ref,choices[start,]$codon,"fs")),
                         by=1:nrow(dtr[values=="-"&length>=2,])]
        deleted=0
        for(i in 1:nrow(deletions)){
            from <- deletions$start[[i]]-deleted
            to   <- deletions$end[[i]]-deleted
            choices<-choices[-c(from:to)]
            deleted <- (deletions$end[[i]] - deletions$start[[i]] +1)
        }
        choices <- rbind(choices,deletions,fill=TRUE)
        setkey(choices,id)
    }

    return(choices)
}

mround <- function(x,base){
    base*round(x/base)
}

report_hetero_indels <- function(calls){
    rev <- !is.null(calls[["call_rev"]])
    if(rev) {
        secondary_seq <- gsub("[ -]","",paste(get_consensus_mut(calls[["mut_call_fwd"]],calls[["mut_call_rev"]],calls[,list(iA_fwd,iC_fwd,iG_fwd,iT_fwd,iA_rev,iC_rev,iG_rev,iT_rev)]),collapse = ""))
    } else secondary_seq <- gsub("[ -]","",paste(calls[["mut_call_fwd"]],collapse = ""))
    primary_seq <- gsub("[ -]","",paste(calls[["user_sample"]],collapse = ""))
    hetero_indel_aln <- pairwiseAlignment(primary_seq, secondary_seq,type = "global-local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
#     # pattern/primary overhang/tail = secondary deletion
#     if(start(pattern(hetero_indel_aln)) > start(subject(hetero_indel_aln))){
#         overhang1 <- paste(replicate(start(pattern(hetero_indel_aln))-1, "+"), collapse = "") # stri_dup("abc",3) from stringi
#     }
#     # subject/secondary overhang/tail = secondary insertion
#     else{
#         overhang1 <- paste(replicate(start(subject(hetero_indel_aln))-1, "-"), collapse = "") # stri_dup("abc",3) from stringi
#     }
#     cmpStr_hetero_indel_aln <- paste(overhang1,compareStrings(hetero_indel_aln),sep="")
#     hetero_del_tab <- stringi::stri_locate_all_regex(cmpStr_hetero_indel_aln,"[\\+]+")[[1]]
#     hetero_ins_tab <- stringi::stri_locate_all_regex(cmpStr_hetero_indel_aln,"[\\-]+")[[1]]
    hetero_ins_tab <- stringi::stri_locate_all_regex(compareStrings(hetero_indel_aln),"[\\-]+")[[1]] + start(pattern(hetero_indel_aln)) - 1
    hetero_del_tab <- stringi::stri_locate_all_regex(compareStrings(hetero_indel_aln),"[\\+]+")[[1]] + start(pattern(hetero_indel_aln)) - 1
    hetero_indel_pid <- round(pid(hetero_indel_aln),1)
    return(list(hetero_indel_aln, hetero_indel_pid, hetero_ins_tab, hetero_del_tab))
}

get_consensus_mut <- function(mut_fwd,mut_rev,intens_tab){
    names <- structure(1:4,names = c("A","C","G","T"))
    pa <- pairwiseAlignment(gsub("[ -]","",paste(mut_fwd,collapse = "")), gsub("[ -]","",paste(mut_rev,collapse = "")),type = "global",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
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

#reconstructs user_mut and mut_peak_pct by shifting user_call_fwd and user_call_rev by detected indels indels
incorporate_hetero_indels_func <- function(calls){
    if(!is.na(g_hetero_del_tab[1])) dels <- as.vector(unlist(apply(g_hetero_del_tab,1,function(x) x[1]:x[2])))
    else dels <- numeric()
    if(!is.na(g_hetero_ins_tab[1])) ins <- as.vector(unlist(apply(g_hetero_ins_tab,1,function(x) x[1]:x[2])))
    else ins <- numeric()
    if(max(length(dels),length(ins)) != 0){
        calls[,het_mut_call_fwd := incorporate_single_vec(calls[["mut_call_fwd"]],ins,dels,"char",T)]
         calls[,het_mut_peak_pct_fwd := incorporate_single_vec(calls[["mut_peak_pct_fwd"]],ins,dels,"num",T)]
#         calls[,het_mut_s2n_abs_fwd := incorporate_single_vec(calls[["mut_s2n_abs_fwd"]],ins,dels,"num",T)]
        calls[set_by_user == FALSE,c("user_mut","mut_peak_pct") := list(het_mut_call_fwd,het_mut_peak_pct_fwd)]
        if(any(colnames(calls) == "call_rev")){
            calls[,het_mut_call_rev := incorporate_single_vec(calls[["mut_call_rev"]],ins,dels,"char",F)]
             calls[,het_mut_peak_pct_rev := incorporate_single_vec(calls[["mut_peak_pct_rev"]],ins,dels,"num",F)]
#             calls[,het_mut_s2n_abs_rev := incorporate_single_vec(calls[["mut_s2n_abs_rev"]],ins,dels,"num",F)]
            calls[(quality_fwd < quality_rev | user_mut == "-") & set_by_user == FALSE,c("user_mut","mut_peak_pct") := list(het_mut_call_rev,het_mut_peak_pct_rev)]
        }
    }
    return(calls)
}

incorporate_single_vec <- function(vec,ins,dels,type,fwd){



    if(type == "num") new_vec <- numeric(length(vec))
    else new_vec <- rep("-",length(vec))

    if(fwd) {
        rem = NULL
        for(i in 1:length(dels)){
            if(g_calls[id==dels[i],mut_call_fwd] =="-")  #picking up a deletion already present in fwd calls
                rem = c(rem,i)
        }
        if(!is.null(rem)) dels <- dels[-rem]
        vec <- vec[1:min(length(vec),length(new_vec) - length(dels))]
        new_vec[setdiff(seq_along(new_vec),dels)] <- vec
    } else {
        vec <- vec[1 + length(vec) -  min(length(vec),length(new_vec) - length(dels)):1]
        new_vec[setdiff(seq_along(new_vec),dels)] <- vec
    }
    if(length(ins) > 0) new_vec <- new_vec[-ins]
    return(new_vec)
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
