library(zoo)
library(seqinr)

# the function extracts the signal intesities for each channel and returns it formatted for the javascript chromatograph
get_intensities <- function(data,calls,deletions,norm=FALSE) {
    #abi file documentation http://www.bioconductor.org/packages/release/bioc/vignettes/sangerseqR/inst/doc/sangerseq_walkthrough.pdf
    intens <- data.table(data$DATA.9,data$DATA.10,data$DATA.11,data$DATA.12)
    
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
    if(nrow(intens)>(calls[length(trace_peak)]$trace_peak+100)){
      intens <- intens[1:(calls[length(trace_peak)]$trace_peak+100)]
    }
    
    intens<-normalize_intensities_lengths(intens,calls[,trace_peak],11)
    fwo <- data$FWO
    intens<-setnames(data.table(intens),c(substring(fwo,1,1),substring(fwo,2,2),substring(fwo,3,3),substring(fwo,4,4)))
    
    #adjust call positions to normalized graph 
    calls <- calls[,trace_peak:=rescale_call_positions(calls[,trace_peak],11)]
    if(length(deletions)!=0){
        intens[,id:=c(1:nrow(intens))]
        setkey(intens,id)
        del_pos <- calls[id %in% deletions][,trace_peak]
        rep = 0
        for(i in c(1:length(del_pos))){
          
                   print(del_pos[i])
                   pos <- del_pos[i]
                   for(i in 1:12){intens<-rbind(intens,list(A=0,C=0,G=0,T=0,id=(pos-6 +rep+ i/100)))} 
                   rep = rep -12
                   #return(intens)
        }
        setkey(intens,id)
        intens[,id:=c(1:nrow(intens))]
        setkey(intens,id)
        
    }
    
    return(list(intens=intens,calls=calls))
}


get_call_data <- function(data,data_r, rm7qual_thres=12, qual_thres=10, aln_min=0.2){
    #TO DO if(length(data$PLOC.1)<=length(data$PBAS.1)){}
    deletions <- list()
    if(is.null(data_r)) {
        qual      <- data$PCON.2
        rm7qual   <- rollmean(qual,k=7)
        
        res       <- generate_ref(data$PBAS.2, aln_min)
        calls <- data.table(id         = seq_along(data$PLOC.2)
                           ,user_mod   = str_split(data$PBAS.2,pattern="")[[1]][seq_along(data$PLOC.2)]
                           ,call       = str_split(data$PBAS.2,pattern="")[[1]][seq_along(data$PLOC.2)]
                           ,reference  = str_split(data$PBAS.2,pattern="")[[1]][seq_along(data$PLOC.2)]
                           ,trace_peak = data$PLOC.2
                           ,quality    = qual)
        calls[,rm7qual := c(quality[1:3],rollmean(quality,k=7),quality[(length(quality) - 2):length(quality)])]
        setkey(res[[2]],id)
        if(length(res[[3]]) > 0) {
            setkey(calls,id)
            add <- calls[as.integer(res[[3]]),]
            data.table::set(add,NULL,"reference",unlist(strsplit(res[[2]][type == "I"][["replace"]],"")))
            calls <- rbind(calls,add[,id := res[[3]]][,call := "-"][,user_mod := "-"])
        }
        deletions <- res[[3]]
        calls <- merge(calls,res[[1]],all.x = T,by = "id")
        data.table::set(calls,which(calls[["id"]] %in% res[[2]][type != "I"][["t_pos"]]),"reference",res[[2]][type != "I"][["replace"]])
        data.table::set(calls,which(is.na(calls[["gen_coord"]])),"reference","NA")
        data.table::set(calls,which(calls[["rm7qual"]] < qual_thres | calls[["quality"]] < qual_thres),"user_mod","low qual")
        
    } else {
        user_align <- get_for_rev_align(data$PBAS.2,data_r$PBAS.2,data$PCON.2,data_r$PCON.2)
        res <- generate_ref(paste(user_align[[1]],collapse = ""), aln_min)
        calls <- data.table(id         = seq_along(user_align[[1]])
                            ,user_mod   = user_align[[1]]
                            ,call       = user_align[[2]]
                            ,call_r     = user_align[[3]]
                            ,reference  =  user_align[[1]]
                            ,trace_peak = data$PLOC.2
                            ,trace_peak_r = data_r$PLOC.2
                            ,quality    = user_align[[4]]
                            ,quality_r    = user_align[[5]])
        cons_qual <- sapply(seq_along(user_align[[4]]),function(x) max(user_align[[4]][x],user_align[[5]][x]))
        calls[,rm7qual := c(cons_qual[1:3],rollmean(cons_qual,k=7),cons_qual[(length(cons_qual) - 2):length(cons_qual)])]
        setkey(res[[2]],id)
        if(length(res[[3]]) > 0) {
            setkey(calls,id)
            add <- calls[as.integer(res[[3]]),]
            data.table::set(add,NULL,"reference",unlist(strsplit(res[[2]][type == "I"][["replace"]],"")))
            calls <- rbind(calls,add[,id := res[[3]]][,call := "-"][,call_r := "-"][,user_mod := "-"])
        }
        calls <- merge(calls,res[[1]],all.x = T,by = "id")
        data.table::set(calls,which(calls[["id"]] %in% res[[2]][type != "I"][["t_pos"]]),"reference",res[[2]][type != "I"][["replace"]])
        data.table::set(calls,which(is.na(calls[["gen_coord"]])),"reference","NA")
        data.table::set(calls,which(calls[["rm7qual"]] < qual_thres | calls[["quality"]] < qual_thres),"user_mod","low qual")
    }
    #helper_intrex contains intesities coordinates of start and end of exons with the sequence id (position in sequence coordinates)
    helperdat <- list()
    helperdat$helper_intrex <- list()
    helperdat$helper_intrex <- setnames(calls[!is.na(exon_intron),list(min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","trace_peak","end"))
    helperdat$helper_intrex <- setnames(merge(helperdat$helper_intrex,calls[,list(id,trace_peak)],by="trace_peak"),"trace_peak","start")
    return(list(calls=calls,helperdat=helperdat,deletions=deletions))
}

generate_ref <-function(user_seq, aln_min=0.2){
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
    for(index in OK_align[order(align$score[OK_align],decreasing = T)]){
        if(any(seq_coverage[(align$start[index] + 1):align$end[index]])) {
            multiple_covered[[ref_names[index]]] <- c(align$start[index],align$end[index],align$score[index])
            OK_align <- OK_align[-which(index == OK_align)]
        } else {
            seq_coverage[(align$start[index] + 1):align$end[index]] <- T
        }
    }

    all_aligns <- seq_along(align$score)
    for(index in all_aligns[order(align$score[all_aligns],decreasing = T)]){
        if(!any(seq_coverage[(align$start[index] + 1):align$end[index]]) &&  align$score[index] / (align$end[index] - align$end[index]) > 0.6 && align$score[index] > 10) {
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

get_for_rev_align <- function(for_seq,rev_seq,for_qual,rev_qual){
    rev_seq <- gsub("[^ACGT]","N",rev_seq)
    for_seq <- gsub("[^ACGT]","N",for_seq)
    sm <- matrix(-1,5,5,dimnames = list(c("A","C","G","T","N"),c("A","C","G","T","N")))
    diag(sm) <- 1
    sm[,"N"] <- 0.1
    sm["N",] <- 0.1
    align <- pairwiseAlignment(pattern = reverseComplement(DNAString(rev_seq)), subject = for_seq,type = "overlap",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
    cons_length <- max(start(subject(align)), start(pattern(align))) - 1 + max(nchar(for_seq) - end(subject(align)), nchar(rev_seq) - end(pattern(align))) + nchar(as.character(subject(align)))
    for_split <- rep("-",cons_length)
    rev_split <- for_split
    splitfr <- strsplit(c(for_seq,as.character(reverseComplement(DNAString(rev_seq)))),"")
    for_split[(1:nchar(as.character(subject(align)))) + max(start(subject(align)), start(pattern(align))) - 1] <- strsplit(as.character(subject(align)),"")[[1]]
    rev_split[(1:nchar(as.character(subject(align)))) + max(start(subject(align)), start(pattern(align))) - 1] <- strsplit(as.character(pattern(align)),"")[[1]]
    if(start(subject(align)) > 1) for_split[1:(start(subject(align)) - 1)] <- splitfr[[1]][1:(start(subject(align)) - 1)]
    if(start(pattern(align)) > 1) rev_split[1:(start(pattern(align)) - 1)] <- splitfr[[2]][1:(start(pattern(align)) - 1)]
    if(nchar(for_seq) - end(subject(align)) > 0) for_split[(length(for_split) - (nchar(for_seq) - end(subject(align))) + 1):length(for_split)] <- splitfr[[1]][length(for_split) - (nchar(for_seq) - end(subject(align))):length(for_split)]
    if(nchar(rev_seq) - end(pattern(align)) > 0) rev_split[(length(rev_split) - (nchar(rev_seq) - end(pattern(align))) + 1):length(rev_split)] <- splitfr[[2]][length(rev_split) - (nchar(rev_seq) - end(pattern(align))):length(rev_split)]
    
    for_split_qual <- rep(0,cons_length)
    rev_split_qual <- for_split_qual
    for_split_qual[which(for_split != "-")] <- for_qual
    rev_split_qual[which(rev_split != "-")] <- rev_qual
    cons_split <- for_split
    cons_split[which(for_split_qual < rev_split_qual)] <- rev_split[which(for_split_qual < rev_split_qual)]
    return(list(cons_split,for_split,rev_split,for_split_qual,rev_split_qual))
}

get_alignment <- function(data,user_seq,cores,type = "overlap"){
    insert_weights <- c(1,4,8,12.6,16.3,25,30,rep(50,100))
    if(length(data) < 1) return(NULL)
    inserts_universal <- F
    if(length(data) <= cores * 2) cores <- 1
    splitid <- sort(seq_along(data) %% cores)
    sm <- matrix(-1,5,5,dimnames = list(c("A","C","G","T","N"),c("A","C","G","T","N")))
    diag(sm) <- 1
    sm[,"N"] <- 0
    sm["N",] <- 0
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
#Adam
normalize_intensities_lengths <- function(intensities, call_positions, intervening_length){
  interpolate <- function(vec, coords, length){
    print(embed(coords,2))
    c(c(vec[1:(coords[1]-1)],
      round(as.vector(apply(embed(coords, 2), 1, function(x) {approx(vec[x[2]:x[1]], n=length + 2)$y[-(length + 2)]})))),
      vec[coords[length(coords)]:length(vec)])
  }
  return(apply(intensities, 2, function(x) interpolate(x, unique(call_positions), intervening_length)))
}
rescale_call_positions <- function(call_positions, intervening_length){
  return(seq(from = call_positions[1], by = intervening_length + 1, length.out = length(call_positions)))
}
#add zero intensities for deletions to be properly shown in graph 
#add_zero_intens(){
  
#}