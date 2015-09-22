library(zoo)           #for the rolling mean function
library(seqinr)

g_base_noise <<- 2

trace_peak_me <- function(p,iA,iC,iG,iT,pA,pC,pG,pT){
    mut_peak <- (sort(c(iA,iC,iG,iT),decreasing = TRUE)[p])
         if (mut_peak == 0 ) return(NA)
    else if (mut_peak == iA) return(pA)
    else if (mut_peak == iC) return(pC)
    else if (mut_peak == iG) return(pG)
    else if (mut_peak == iT) return(pT)
}

trace_peak_me2 <- function(i,p){
    return(p[which.max(i)])
}

get_call_data <- function(data, data_rev){
    #TO DO if(length(data$PLOC.1)<=length(data$PBAS.1)){}
    deletions <- list()
    if(is.null(data_rev)) {
        pri   <- as.character(data@primarySeq)
        res   <- generate_ref(pri)
        calls <- data.table(id        = seq_along(data@primarySeq)
                           ,user_sample  = str_split(data@primarySeq,pattern="")[[1]][seq_along(data@primarySeq)]
                           ,call      = str_split(data@primarySeq,pattern="")[[1]][seq_along(data@primarySeq)]
                           ,reference = str_split(data@primarySeq,pattern="")[[1]][seq_along(data@primarySeq)]
                           ,pA_fwd    = data@peakPosMatrix[,1],pC_fwd = data@peakPosMatrix[,2],pG_fwd = data@peakPosMatrix[,3],pT_fwd = data@peakPosMatrix[,4]
                           ,iA_fwd    = data@peakAmpMatrix[,1],iC_fwd = data@peakAmpMatrix[,2],iG_fwd = data@peakAmpMatrix[,3],iT_fwd = data@peakAmpMatrix[,4]
        )
        intens_calls <- get_intens_calls(pri,data@peakAmpMatrix)
        calls[,trace_peak := trace_peak_me(1,iA_fwd,iC_fwd,iG_fwd,iT_fwd,pA_fwd,pC_fwd,pG_fwd,pT_fwd),by=1:nrow(calls)]
        setkey(res[[2]],id)
        if(length(res[[3]]) > 0) {
            setkey(calls,id)
            add <- calls[as.integer(res[[3]]),]
            data.table::set(add,NULL,"reference",unlist(strsplit(res[[2]][type == "I"][["replace"]],"")))
            calls <- rbind(calls,add[,id := res[[3]]][,call := "-"][,user_sample := "-"])
        }
        trace_peak_fwd <- sapply(1:nrow(data@peakAmpMatrix),function(x) trace_peak_me2(data@peakAmpMatrix[x,],data@peakPosMatrix[x,]))
        trace_peak_rev <- NULL

    } else {

        user_align <- get_fwd_rev_align(as.character(data@primarySeq),as.character(data_rev@primarySeq))
        intens_calls <- get_intens_calls(user_align[[1]],data@peakAmpMatrix,user_align[[2]],data_rev@peakAmpMatrix)
        consensus_seq <- create_consensus_seq(user_align[[1]],user_align[[2]],intens_calls)
        res <- generate_ref(paste(consensus_seq,collapse = ""))

        calls <- data.table(id         = seq_along(consensus_seq)
                            ,user_sample  = consensus_seq
                            ,call      = user_align[[1]]
                            ,call_rev  = user_align[[2]]
                            ,reference = consensus_seq
        )
        calls <- cbind(calls,intens_calls)
        trace_peak_fwd <- sapply(1:nrow(data@peakAmpMatrix),function(x) trace_peak_me2(data@peakAmpMatrix[x,],data@peakPosMatrix[x,]))
        trace_peak_rev <- sapply(1:nrow(data_rev@peakAmpMatrix),function(x) trace_peak_me2(data_rev@peakAmpMatrix[x,],data_rev@peakPosMatrix[x,]))
        setkey(res[[2]],id)
        if(length(res[[3]]) > 0) {
            setkey(calls,id)
            add <- calls[as.integer(res[[3]]),]
            data.table::set(add,NULL,"reference",unlist(strsplit(res[[2]][type == "I"][["replace"]],"")))
            calls <- rbind(calls,add[,id := res[[3]]][,call := "-"][,call_rev := "-"][,user_sample := "-"])
        }
    }

    deletions <- res[[3]]
    calls <- merge(calls,res[[1]],all.x = T,by = "id")

    data.table::set(calls,which(calls[["id"]] %in% res[[2]][type != "I"][["t_pos"]]),"reference",res[[2]][type != "I"][["replace"]])
    data.table::set(calls,which(is.na(calls[["gen_coord"]])),"reference","NA")

    return(list(calls=calls,deletions=deletions,trace_peak_fwd = trace_peak_fwd,trace_peak_rev = trace_peak_rev))
}

create_consensus_seq <- function(fwd_seq,rev_seq,intens_tab){
    fwd_noise <- apply(intens_tab[,1:4,with = F],1,get_noise)
    rev_noise <- apply(intens_tab[,5:8,with = F],1,get_noise)
    consensus_seq <- fwd_seq
    consensus_seq[which(fwd_noise > rev_noise)] <- rev_seq[which(fwd_noise > rev_noise)]
    return(consensus_seq)
}

get_noise <- function(row){
    return((mean(row[order(row)[c(1,2)]]) + g_base_noise) / (row[which.max(row)]  + g_base_noise))
}

get_intens_calls <- function(calls_fwd,intens_fwd,calls_rev=NULL,intens_rev=NULL){
    calls_intens_fwd <- numeric(length(calls_fwd))
    res_fwd <- lapply(1:4,function(x) {
        calls_intens_fwd[which(calls_fwd != "-")] <- intens_fwd[,x]
        return(calls_intens_fwd)
        })
    if(!is.null(calls_rev)){
        calls_intens_rev <- numeric(length(calls_rev))
        res_rev <- lapply(4:1,function(x) {
            calls_intens_rev[which(calls_rev != "-")] <- rev(intens_rev[,x])
            return(calls_intens_rev)
        })
        return(data.table(iA_fwd = res_fwd[[1]],iC_fwd = res_fwd[[2]],iG_fwd = res_fwd[[3]],iT_fwd = res_fwd[[4]],iA_rev = res_rev[[1]],iC_rev = res_rev[[2]],iG_rev = res_rev[[3]],iT_rev = res_rev[[4]]))

    }
    return(data.table(iA_fwd = res_fwd[[1]],iC_fwd = res_fwd[[2]],iG_fwd = res_fwd[[3]],iT_fwd = res_fwd[[4]]))

}

generate_ref <-function(user_seq){
    cores <- 2
    multiple_covered <- list()
    user_seq <- gsub("[^ACGT]","N",user_seq)

    refs   <- readLines("../../data/ref_ex_in.fa")
    ref_info <- gsub(">ref_","",perl = T,refs[seq(1,length(refs),2)])
    ref_info <- strsplit(ref_info,split = "_")
    ref_names <- sapply(ref_info,function(x) x[1])
    ref_start <- as.numeric(sapply(ref_info,function(x) x[2]))
    ref_end <- as.numeric(sapply(ref_info,function(x) x[3]))
    refs <- DNAStringSet(refs[seq(2,length(refs),2)])
    align <- get_alignment(refs,user_seq,cores)
    OK_align <- which(align$score / nchar(refs) > 0.8 | align$score > 50)

    seq_coverage <- logical(nchar(user_seq))
    for(index in OK_align[order(-(align$ref_starts + nchar(refs) - align$ref_ends - 1)[OK_align],align$score[OK_align],decreasing = T)]){
        if(any(seq_coverage[(align$start[index] + 1):align$end[index]])) {
            indeces <- which(!seq_coverage[(align$start[index] + 1):align$end[index]]) + align$start[index]
            if(length(indeces > 50)) {
                align$starts[index] <- min(indeces) - 1
                align$ends[index] <- max(indeces)
                ref_indeces <- which(!seq_coverage[(align$start[index] + 1):align$end[index]]) + align$ref_start[index] - 1
                align$ref_start[index] <- min(ref_indeces)
                align$ref_end[index] <- max(ref_indeces)
                seq_coverage[indeces] <- T
            } else {
                OK_align <- OK_align[-which(index == OK_align)]
            }
            multiple_covered[[ref_names[index]]] <- c(align$start[index],align$end[index],align$score[index])
        } else {
            seq_coverage[(align$start[index] + 1):align$end[index]] <- T
        }
    }

    all_aligns <- seq_along(align$score)
    for(index in all_aligns[order(align$score[all_aligns],decreasing = T)]){
        if(!any(seq_coverage[(align$start[index] + 1):align$end[index]]) && align$score[index] / (abs(align$end[index] - align$start[index])) > 0.6 && align$score[index] > 10) {
            seq_coverage[(align$start[index] + 1):align$end[index]] <- T
            OK_align <- c(OK_align,index)
        }
    }

    user_seq_vs_genome <- rbindlist(lapply(seq_along(OK_align),function(x) data.table(exon_intron = ref_names[OK_align][x]
		,id = sort(c((align$start[OK_align][x] + 1):align$end[OK_align][x],add_insert(align$diffs[type == "I" & id == OK_align[x]])))
		,gen_coord = get_coord(align$start[OK_align][x],align$ref_start[OK_align][x],align$ref_end[OK_align][x],ref_start[OK_align][x],ref_end[OK_align][x],align$diffs[type == "D" & id == OK_align[x]]))))

    align$diffs[type == "D",t_pos := t_pos + 1]
    align$diffs <- align$diffs[id %in% OK_align,]
    return(list(user_seq_vs_genome,align$diffs,add_insert(align$diffs[type == "I"])))
}

add_insert <- function(diffs){
    if(nrow(diffs) > 0) return(unlist(lapply(1:nrow(diffs),function(x) diffs[["t_pos"]][x] + (1:nchar(diffs[["replace"]][x]))/10 )))
    else return(numeric())
}

get_coord <- function(seq_start,al_start,al_end,ref_start,ref_end,diffs){
    coord <- (ref_start:ref_end)[al_start:al_end]
    if(nrow(diffs) > 0){
        diffs[,t_pos := t_pos - seq_start + 1]
        coord_out <- numeric(length(coord) + nrow(diffs))
        coord_out[-diffs[["t_pos"]]] <- coord
        for(index in which(coord_out == 0)) coord_out[index] <- coord_out[index - 1] + 0.1
        return(coord_out)
    } else return(coord)
}

get_fwd_rev_align <- function(fwd_seq,rev_seq){
    #rev_seq <- gsub("[^ACGT]","N",rev_seq)
    #fwd_seq <- gsub("[^ACGT]","N",fwd_seq)
    #sm <- matrix(-1,5,5,dimnames = list(c("A","C","G","T","N"),c("A","C","G","T","N")))
    #diag(sm) <- 1
    #sm[,"N"] <- 0.1
    #sm["N",] <- 0.1
    #sm <- matrix(c(5 ,-4 ,-4 ,-4 ,-4 ,1 ,1 ,-4 ,-4 ,1 ,-4 ,-1 ,-1 ,-1 ,-2 ,-4 ,5 ,-4 ,-4 ,-4 ,1 ,-4 ,1 ,1 ,-4 ,-1 ,-4 ,-1 ,-1 ,-2 ,-4 ,-4 ,5 ,-4 ,1 ,-4 ,1 ,-4 ,1 ,-4 ,-1 ,-1 ,-4 ,-1 ,-2 ,-4 ,-4 ,-4 ,5 ,1 ,-4 ,-4 ,1 ,-4 ,1 ,-1 ,-1 ,-1 ,-4 ,-2 ,-4 ,-4 ,1 ,1 ,-1 ,-4 ,-2 ,-2 ,-2 ,-2 ,-1 ,-1 ,-3 ,-3 ,-1 ,1 ,1 ,-4 ,-4 ,-4 ,-1 ,-2 ,-2 ,-2 ,-2 ,-3 ,-3 ,-1 ,-1 ,-1 ,1 ,-4 ,1 ,-4 ,-2 ,-2 ,-1 ,-4 ,-2 ,-2 ,-3 ,-1 ,-3 ,-1 ,-1 ,-4 ,1 ,-4 ,1 ,-2 ,-2 ,-4 ,-1 ,-2 ,-2 ,-1 ,-3 ,-1 ,-3 ,-1 ,-4 ,1 ,1 ,-4 ,-2 ,-2 ,-2 ,-2 ,-1 ,-4 ,-1 ,-3 ,-3 ,-1 ,-1 ,1 ,-4 ,-4 ,1 ,-2 ,-2 ,-2 ,-2 ,-4 ,-1 ,-3 ,-1 ,-1 ,-3 ,-1 ,-4 ,-1 ,-1 ,-1 ,-1 ,-3 ,-3 ,-1 ,-1 ,-3 ,-1 ,-2 ,-2 ,-2 ,-1 ,-1 ,-4 ,-1 ,-1 ,-1 ,-3 ,-1 ,-3 ,-3 ,-1 ,-2 ,-1 ,-2 ,-2 ,-1 ,-1 ,-1 ,-4 ,-1 ,-3 ,-1 ,-3 ,-1 ,-3 ,-1 ,-2 ,-2 ,-1 ,-2 ,-1 ,-1 ,-1 ,-1 ,-4 ,-3 ,-1 ,-1 ,-3 ,-1 ,-3 ,-2 ,-2 ,-2 ,-1 ,-1 ,-2 ,-2 ,-2 ,-2 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1),15,15,dimnames = list(c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N")))
    sm <- matrix(c(1 ,-1 ,-1 ,-1 ,-1 ,0.5 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,0.1 ,0 ,-1 ,1 ,-1 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.5 ,-1 ,0.1 ,-1 ,0.1 ,0.1 ,0 ,-1 ,-1 ,1 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,-1 ,0.1 ,0 ,-1 ,-1 ,-1 ,1 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.1 ,0.1 ,0.1 ,-1 ,0 ,-1 ,-1 ,0.5 ,0.5 ,0.1 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.5 ,0.5 ,-1 ,-1 ,-1 ,0.1 ,0 ,0 ,0 ,0 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.1 ,0.5 ,-1 ,0.5 ,-1 ,0 ,0 ,0.1 ,-1 ,0 ,0 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-1 ,0.5 ,-1 ,0.5 ,0 ,0 ,-1 ,0.1 ,0 ,0 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.5 ,0.5 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,-1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.5 ,-1 ,-1 ,0.5 ,0 ,0 ,0 ,0 ,-1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0 ,0.1 ,0 ,0 ,0.1 ,0.1 ,0.1 ,-1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0.1 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,-1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1),15,15,dimnames = list(c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N")))

    align <- pairwiseAlignment(pattern = reverseComplement(DNAString(rev_seq)), subject = fwd_seq,type = "overlap",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
    cons_length <- max(start(subject(align)), start(pattern(align))) - 1 + max(nchar(fwd_seq) - end(subject(align)), nchar(rev_seq) - end(pattern(align))) + nchar(as.character(subject(align)))
    fwd_split <- rep("-",cons_length)
    rev_split <- fwd_split
    splitfr <- strsplit(c(fwd_seq,as.character(reverseComplement(DNAString(rev_seq)))),"")
    fwd_split[(1:nchar(as.character(subject(align)))) + max(start(subject(align)), start(pattern(align))) - 1] <- strsplit(as.character(subject(align)),"")[[1]]
    rev_split[(1:nchar(as.character(subject(align)))) + max(start(subject(align)), start(pattern(align))) - 1] <- strsplit(as.character(pattern(align)),"")[[1]]
    if(start(subject(align)) > 1) fwd_split[1:(start(subject(align)) - 1)] <- splitfr[[1]][1:(start(subject(align)) - 1)]
    if(start(pattern(align)) > 1) rev_split[1:(start(pattern(align)) - 1)] <- splitfr[[2]][1:(start(pattern(align)) - 1)]
    if(nchar(fwd_seq) - end(subject(align)) > 0) fwd_split[(length(fwd_split) - (nchar(fwd_seq) - end(subject(align))) + 1):length(fwd_split)] <- splitfr[[1]][(end(subject(align)) + 1):length(splitfr[[1]])]
    if(nchar(rev_seq) - end(pattern(align)) > 0) rev_split[(length(rev_split) - (nchar(rev_seq) - end(pattern(align))) + 1):length(rev_split)] <- splitfr[[2]][(end(pattern(align)) + 1):length(splitfr[[2]])]

    fwd_offset <- start(pattern(align))
    rev_offset <- start(subject(align))



    return(list(fwd_split,rev_split,fwd_offset,rev_offset))
}

get_alignment <- function(data,user_seq,cores,type = "overlap"){
    insert_weights <- c(1,4,8,12.6,16.3,25,30,rep(50,100))
    if(length(data) < 1) return(NULL)
    inserts_universal <- F
    if(length(data) <= cores * 2) cores <- 1
    splitid <- sort(seq_along(data) %% cores)
    # <- matrix(-1,5,5,dimnames = list(c("A","C","G","T","N"),c("A","C","G","T","N")))
    #diag(sm) <- 1
    #sm[,"N"] <- 0
    #sm["N",] <- 0
    #sm <- matrix(c(5 ,-4 ,-4 ,-4 ,-4 ,1 ,1 ,-4 ,-4 ,1 ,-4 ,-1 ,-1 ,-1 ,-2 ,-4 ,5 ,-4 ,-4 ,-4 ,1 ,-4 ,1 ,1 ,-4 ,-1 ,-4 ,-1 ,-1 ,-2 ,-4 ,-4 ,5 ,-4 ,1 ,-4 ,1 ,-4 ,1 ,-4 ,-1 ,-1 ,-4 ,-1 ,-2 ,-4 ,-4 ,-4 ,5 ,1 ,-4 ,-4 ,1 ,-4 ,1 ,-1 ,-1 ,-1 ,-4 ,-2 ,-4 ,-4 ,1 ,1 ,-1 ,-4 ,-2 ,-2 ,-2 ,-2 ,-1 ,-1 ,-3 ,-3 ,-1 ,1 ,1 ,-4 ,-4 ,-4 ,-1 ,-2 ,-2 ,-2 ,-2 ,-3 ,-3 ,-1 ,-1 ,-1 ,1 ,-4 ,1 ,-4 ,-2 ,-2 ,-1 ,-4 ,-2 ,-2 ,-3 ,-1 ,-3 ,-1 ,-1 ,-4 ,1 ,-4 ,1 ,-2 ,-2 ,-4 ,-1 ,-2 ,-2 ,-1 ,-3 ,-1 ,-3 ,-1 ,-4 ,1 ,1 ,-4 ,-2 ,-2 ,-2 ,-2 ,-1 ,-4 ,-1 ,-3 ,-3 ,-1 ,-1 ,1 ,-4 ,-4 ,1 ,-2 ,-2 ,-2 ,-2 ,-4 ,-1 ,-3 ,-1 ,-1 ,-3 ,-1 ,-4 ,-1 ,-1 ,-1 ,-1 ,-3 ,-3 ,-1 ,-1 ,-3 ,-1 ,-2 ,-2 ,-2 ,-1 ,-1 ,-4 ,-1 ,-1 ,-1 ,-3 ,-1 ,-3 ,-3 ,-1 ,-2 ,-1 ,-2 ,-2 ,-1 ,-1 ,-1 ,-4 ,-1 ,-3 ,-1 ,-3 ,-1 ,-3 ,-1 ,-2 ,-2 ,-1 ,-2 ,-1 ,-1 ,-1 ,-1 ,-4 ,-3 ,-1 ,-1 ,-3 ,-1 ,-3 ,-2 ,-2 ,-2 ,-1 ,-1 ,-2 ,-2 ,-2 ,-2 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1 ,-1),15,15,dimnames = list(c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N")))
    sm <- matrix(c(1 ,-1 ,-1 ,-1 ,-1 ,0.5 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,0.1 ,0 ,-1 ,1 ,-1 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.5 ,-1 ,0.1 ,-1 ,0.1 ,0.1 ,0 ,-1 ,-1 ,1 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,-1 ,0.1 ,0 ,-1 ,-1 ,-1 ,1 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.1 ,0.1 ,0.1 ,-1 ,0 ,-1 ,-1 ,0.5 ,0.5 ,0.1 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.5 ,0.5 ,-1 ,-1 ,-1 ,0.1 ,0 ,0 ,0 ,0 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.1 ,0.5 ,-1 ,0.5 ,-1 ,0 ,0 ,0.1 ,-1 ,0 ,0 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-1 ,0.5 ,-1 ,0.5 ,0 ,0 ,-1 ,0.1 ,0 ,0 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.5 ,0.5 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,-1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.5 ,-1 ,-1 ,0.5 ,0 ,0 ,0 ,0 ,-1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0 ,0.1 ,0 ,0 ,0.1 ,0.1 ,0.1 ,-1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0.1 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,-1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1),15,15,dimnames = list(c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N")))

    a <- lapply(1:cores - 1L, function(x) {
        align <- pairwiseAlignment(pattern = data[which(splitid == x)], subject = user_seq,type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -0.3)
        res <- list()
        res$starts <- start(subject(align)) - 1
        res$ends <- end(subject(align))
        res$ref_starts <- start(pattern(align))
        res$ref_ends <- end(pattern(align))
        mismatches <- as.data.table(mismatchTable(align))
        setnames(mismatches,c("PatternId","SubjectStart","PatternSubstring","SubjectSubstring"),c("id","t_pos","replace","orig"))
        mismatches[,c("PatternEnd","PatternStart","SubjectEnd") := NULL]
        mismatches[,type := "M"]
        setcolorder(mismatches, c(1,5,3,4,2))
        mismatches[,weight := 1]
        insert_ids <- sapply(insertion(align),length)
        if(sum(insert_ids) > 0){
            insert_start <- unlist(lapply(insertion(align)[which(insert_ids > 0)],start))
            insert_width <- unlist(lapply(insertion(align)[which(insert_ids > 0)],width))
            if(!inserts_universal){
                insert_offset <- unlist(lapply(insertion(align)[which(insert_ids > 0)],function(x) cumsum(c(0,width(x)[-length(width(x))]))))
                insert_replace <- stri_sub(as.character(rep(aligned(pattern(align)),insert_ids)),insert_start + insert_offset,length = insert_width)
                insertions <- data.table(id = rep(1:length(insert_ids),insert_ids),type = "I",t_pos = insert_start + rep(res$starts,insert_ids),orig = "-",replace = insert_replace,weight = insert_weights[insert_width])
            } else {
                insertions <- data.table(id = rep(1:length(insert_ids),insert_ids),type = "I",t_pos = insert_start + rep(res$starts,insert_ids),orig = "-",replace = as.character(insert_width),weight = insert_weights[insert_width])
            }
        } else insertions <- NULL
        del_ids <- sapply(deletion(align),function(x) sum(width(x)))
        if(sum(del_ids) > 0){
            del_pos <- unlist(gregexpr("\\-",gsub("\\+","",compareStrings(align))))
            split_t <- strsplit(as.character(user_seq),"")[[1]]
            deletions <- data.table(id = rep(1:length(del_ids),del_ids),type = "D",t_pos = del_pos[del_pos>0] + rep(res$starts,del_ids),orig = split_t[del_pos[del_pos>0] + rep(res$starts,del_ids)],replace = "-",weight = 1)
        } else deletions <- NULL
        res$diffs <- rbindlist(list(mismatches,insertions,deletions))
        res$diffs[,id := id + length(which(splitid < x))]
        res$score <- Biostrings::score(align)
        return(res)
    })

    #,mc.preschedule = F,mc.cores = cores

    res <- list()
    res$score <- unlist(sapply(a,function(x) x$score))
    res$ref_starts <- unlist(sapply(a,function(x) x$ref_starts))
    res$ref_ends <- unlist(sapply(a,function(x) x$ref_ends))
    res$starts <- unlist(sapply(a,function(x) x$starts))
    res$ends <- unlist(sapply(a,function(x) x$ends))
    res$diffs <- rbindlist(lapply(a,function(x) x$diffs))
    res$diffs[,orig := as.character(orig)]
    res$diffs[,replace := as.character(replace)]
    set(res$diffs,which(res$diffs[["replace"]] == "1"),"replace","I")
    return(res)
}

# the function extracts the signal intesities for each channel
get_intensities <- function(data,data_rev=NULL,trace_peak_fwd,trace_peak_rev,calls,deletions=NULL,norm=FALSE) {
    #abi file documentation http://www.bioconductor.org/packages/release/bioc/vignettes/sangerseqR/inst/doc/sangerseq_walkthrough.pdf
    rev <- !(is.null(data_rev))

            intens     <- data.table(data@traceMatrix[,1],data@traceMatrix[,2],data@traceMatrix[,3],data@traceMatrix[,4])
    if(rev) intens_rev <- data.table(data_rev@traceMatrix[,1],data_rev@traceMatrix[,2],data_rev@traceMatrix[,3],data_rev@traceMatrix[,4])
    else    intens_rev <- NULL

    #!!out of date
    if(norm==TRUE){
        #first and last 500 points
        f_intens_start <- data.table()
        f_intens_end   <- data.table()
        last <- nrow(intens)
        for(i in 1:499){
            f_intens_start <- rbind(f_intens_start,as.list(apply(intens[1:(i+499)],2,function (x){sum(x)/(499+i)})))
            f_intens_end   <- rbind(f_intens_end,as.list(apply(intens[(last-998+i):last],2,function(x){sum(x)/(998-i)})))
        }
        #1000 rolling mean
        f_intens_mid <- data.table(rollmean(intens,k=999))
        f_intens <- rbind(f_intens_start,f_intens_mid,f_intens_end)
        intens<-intens/f_intens
    }

    #cliping the end of chromatogram after last call
    if(nrow(intens) > max(trace_peak_fwd)){
        intens <- intens[1:max(trace_peak_fwd)]
    }
    #cliping the end of reverse chromatogram after last call + length of first trace peak in fwd so we get the same offset once we turn it over
    offset <- trace_peak_fwd[1]
    if(rev){
        if(nrow(intens_rev) > max(trace_peak_rev)+offset){
            intens_rev <- intens_rev[1:max(trace_peak_rev)+offset]
        }
    }


    intens <- normalize_peak_width(intens,trace_peak_fwd,11)
    intens <- setnames(data.table(intens),c("A","C","G","T"))
    calls  <- calls[,trace_peak:=rescale_call_positions(trace_peak_fwd[1],nrow(calls),11)]
    if(rev){
        intens_rev <- normalize_peak_width(intens_rev,trace_peak_rev,11)
        intens_rev <- setnames(data.table(intens_rev),c("T","G","C","A"))
        intens_rev <- intens_rev[nrow(intens_rev):1]
        calls      <- calls[,trace_peak_rev:=rescale_call_positions(trace_peak_rev[1],nrow(calls),11)]
    }

    if(rev){
        deletions     <- calls[call    =="-"][,id]
        deletions_rev <- calls[call_rev=="-"][,id]
    }else deletions_rev <- list()
    if(length(deletions)!=0){
        intens[,id:=c(1:nrow(intens))]
        setkey(intens,id)
        del_pos <- 0
        del_pos <- calls[id %in% deletions][,trace_peak]
        rep <- 0
        for(i in c(1:length(del_pos))){
            pos <- del_pos[i]
            for(i in 1:12){intens<-rbind(intens,list(A=0,C=0,G=0,T=0,id=(pos-6 +rep+ i/100)))}
            rep <- rep -12
        }
        setkey(intens,id)
        intens[,id:=c(1:nrow(intens))]
        setkey(intens,id)
    }
    if(length(deletions_rev)!=0){
        intens_rev[,id:=c(1:nrow(intens_rev))]
        setkey(intens_rev,id)
        del_pos <- 0
        del_pos <- calls[id %in% deletions_rev][,trace_peak_rev]
        rep <- 0
        for(i in c(1:length(del_pos))){
            pos <- del_pos[i]
            for(i in 1:12){intens_rev<-rbind(intens_rev,list(A=0,C=0,G=0,T=0,id=(pos-6 +rep+ i/100)))}
            rep <- rep -12
        }
        setkey(intens_rev,id)
        intens_rev[,id:=c(1:nrow(intens_rev))]
        setkey(intens_rev,id)
    }

    return(list(intens=intens,intens_rev=intens_rev,calls=calls))
}

#Adam
normalize_peak_width <- function(intensities, call_positions, intervening_length){
    interpolate <- function(vec, coords, length){
        c(c(vec[1:(coords[1]-1)],
        round(as.vector(apply(embed(coords, 2), 1, function(x) {approx(vec[x[2]:x[1]], n=length + 2)$y[-(length + 2)]})))),
        vec[coords[length(coords)]:length(vec)])
    }
    return(apply(intensities, 2, function(x) interpolate(x, unique(call_positions), intervening_length)))
}
rescale_call_positions <- function(call_positions_start, call_positions_length, intervening_length){
    return(seq(from = call_positions_start, by = intervening_length + 1, length.out = call_positions_length))
}
