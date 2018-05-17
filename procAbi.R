library(zoo)       # for the rolling mean function
library(seqinr)    # https://cran.r-project.org/web/packages/seqinr/seqinr.pdf

# g_base_noise <<- 2
# g_calibration_length <<- 30

sm <<- matrix(c(1 ,-1 ,-1 ,-1 ,-1 ,0.5 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,0.1 ,0 ,-1 ,1 ,-1 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.5 ,-1 ,0.1 ,-1 ,0.1 ,0.1 ,0 ,-1 ,-1 ,1 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,-1 ,0.1 ,0 ,-1 ,-1 ,-1 ,1 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.1 ,0.1 ,0.1 ,-1 ,0 ,-1 ,-1 ,0.5 ,0.5 ,0.1 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.5 ,0.5 ,-1 ,-1 ,-1 ,0.1 ,0 ,0 ,0 ,0 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.1 ,0.5 ,-1 ,0.5 ,-1 ,0 ,0 ,0.1 ,-1 ,0 ,0 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-1 ,0.5 ,-1 ,0.5 ,0 ,0 ,-1 ,0.1 ,0 ,0 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.5 ,0.5 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,-1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.5 ,-1 ,-1 ,0.5 ,0 ,0 ,0 ,0 ,-1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0 ,0.1 ,0 ,0 ,0.1 ,0.1 ,0.1 ,-1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0.1 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,-1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1),15,15,dimnames = list(c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N")))

get_call_data <- function(data,data_rev,single_rev,glassed_ref){
    #TO DO if(length(data$PLOC.1)<=length(data$PBAS.1)){}
    deletions <- list()
    if(is.null(data_rev)) {
        qual  <- data$PCON.1
        if(is.null(qual)){
            qual <- get_pseudo_qual(data)
            qual_present <- FALSE
        }else{
            qual_present <- TRUE
        }
        if(single_rev) {
            res   <- generate_ref(complement(reverse(gsub('[\\*,!]','',data$PBAS.1))),glassed_ref)
            calls <- data.table(id           = seq_along(data$PLOC.1)
                                ,user_sample = str_split(complement(gsub('[\\*,!]','',data$PBAS.1)),pattern="")[[1]][seq_along(data$PLOC.1)][length(seq_along(data$PLOC.1)):1]
                                ,call        = str_split(complement(gsub('[\\*,!]','',data$PBAS.1)),pattern="")[[1]][seq_along(data$PLOC.1)][length(seq_along(data$PLOC.1)):1]
                                ,reference   = str_split(complement(gsub('[\\*,!]','',data$PBAS.1)),pattern="")[[1]][seq_along(data$PLOC.1)][length(seq_along(data$PLOC.1)):1]
                                ,trace_peak  = data$PLOC.1
                                ,quality     = qual[length(qual):1]
                                ,VAF         = 0)
        }else{
            res   <- generate_ref(gsub('[\\*,!]','',toupper(substr(data$PBAS.1,0,length(data$PLOC.1)))),glassed_ref)
            calls <- data.table(id           = seq_along(data$PLOC.1)
                                ,user_sample = str_split(gsub('[\\*,!]','',substr(data$PBAS.1,0,length(data$PLOC.1))),pattern="")[[1]][seq_along(data$PLOC.1)]
                                ,call        = str_split(gsub('[\\*,!]','',substr(data$PBAS.1,0,length(data$PLOC.1))),pattern="")[[1]][seq_along(data$PLOC.1)]
                                ,reference   = str_split(gsub('[\\*,!]','',toupper(substr(data$PBAS.1,0,length(data$PLOC.1)))),pattern="")[[1]][seq_along(data$PLOC.1)]
                                ,trace_peak  = data$PLOC.1
                                ,quality     = qual
                                ,quality_fwd = qual
                                ,VAF         = 0)
        }
        # calls[,rm7qual := c(quality[1:3],rollmean(quality,k=7),quality[(length(quality) - 2):length(quality)])]
        setkey(res[[2]],t_pos)
        if(length(res[[3]]) > 0) {
            setkey(calls,id)
            add <- calls[as.integer(res[[3]]),]
            data.table::set(add,NULL,"reference",unlist(strsplit(res[[2]][type == "I"][["replace"]],"")))
            calls <- rbind(calls,add[,id := res[[3]]][,call := "-"][,user_sample := "-"])
        }

    } else {
        fwd_qual <- data$PCON.2
        rev_qual <- data_rev$PCON.2
        qual_present <- TRUE
        if(is.null(fwd_qual)){
            fwd_qual <- get_pseudo_qual(data)
            qual_present <- FALSE
        }
        if(is.null(rev_qual)){
            rev_qual <- get_pseudo_qual(rev_data)
            qual_present <- FALSE
        }
        user_align <- get_fwd_rev_align(gsub('[\\*,!]','',data$PBAS.1),gsub('[\\*,!]','',data_rev$PBAS.1),fwd_qual,rev_qual)
        res <- generate_ref(paste(user_align[[1]],collapse = ""),glassed_ref)
        calls <- data.table(id              = seq_along(user_align[[1]])
                            ,user_sample    = user_align[[1]]
                            ,call           = user_align[[2]]
                            ,call_rev       = user_align[[3]]
                            ,reference      = user_align[[1]]
                            ,quality        = sapply(seq_along(user_align[[4]]),function(x) max(user_align[[4]][x],user_align[[5]][x]))
                            ,quality_fwd    = user_align[[4]]
                            ,quality_rev    = user_align[[5]]
                            ,trace_peak     = data$PLOC.1
                            ,trace_peak_rev = data_rev$PLOC.1
                            ,VAF            = 0)
#         calls[,rm7qual := c(quality[1:3],rollmean(quality,k=7),quality[(length(quality) - 2):length(quality)])]
#         calls[,rm7qual_fwd := c(quality_fwd[1:3],rollmean(quality_fwd,k=7),quality[(length(quality_fwd) - 2):length(quality_fwd)])]
#         calls[,rm7qual_rev := c(quality_rev[1:3],rollmean(quality_rev,k=7),quality[(length(quality_rev) - 2):length(quality_rev)])]
        ##setkey(res[[2]],id) ### sample 2168 id is duplicite??
        setkey(res[[2]],t_pos) ### CHECK!!!
        if(length(res[[3]]) > 0) {
            setkey(calls,id)
            add <- calls[as.integer(res[[3]]),]
            data.table::set(add,NULL,"reference",unlist(strsplit(res[[2]][type == "I"][["replace"]],"")))
            calls <- rbind(calls,add[,id := res[[3]]][,call := "-"][,call_rev := "-"][,user_sample := "-"])
        }
    }

    deletions <- res[[3]]
    calls <- merge(calls,res[[1]],all.x = T,by = "id")
    if(nrow(res[[2]])>0){
        data.table::set(calls,which(calls[["id"]] %in% res[[2]][type != "I"][["t_pos"]]),"reference",res[[2]][type != "I"][["replace"]])
    }
    if(!gsub("data/refs/","",gsub(".glassed.intrex.fasta","",glassed_ref))=="-"){
        data.table::set(calls,which(is.na(calls[["gen_coord"]])),"reference","NA")
    }
    if(! FALSE %in% is.na(calls$exon_intron )) {
        calls[,exon_intron := NULL]
        calls[,exon_intron := 'none']
    }
    if(nrow(calls)<=11){
        stop("none or too few called bases in abi file.")
    }
    return(list(calls=calls,deletions=deletions,qual_present=qual_present))
}

# create_consensus_seq <- function(fwd_seq,rev_seq,intens_tab){
#     intens_tab <- copy(intens_tab)
#     fwd_noise <- apply(intens_tab[,1:4,with = F],1,get_noise)
#     rev_noise <- apply(intens_tab[,5:8,with = F],1,get_noise)
#     fwd_noise[min(which(fwd_noise < 1)):(min(which(fwd_noise < 1)) + g_calibration_length)] <- 1
#     rev_noise[max(which(rev_noise < 1)):(max(which(rev_noise < 1)) - g_calibration_length)] <- 1
#     consensus_seq <- fwd_seq
#     consensus_seq[which(fwd_noise > rev_noise)] <- rev_seq[which(fwd_noise > rev_noise)]
#     return(consensus_seq)
# }

# get_noise <- function(row){
#     return((mean(row[-which.max(row)]) + g_base_noise) / (row[which.max(row)]  + g_base_noise))
# }


get_pseudo_qual <- function(data){
    to = length(data$PLOC.1)
    qual <- c(
        round((1-noise(max(data$DATA.9[c((data$PLOC.1[1]-3):(data$PLOC.1[1]+3))]),
                       max(data$DATA.10[c((data$PLOC.1[1]-3):(data$PLOC.1[1]+3))]),
                       max(data$DATA.11[c((data$PLOC.1[1]-3):(data$PLOC.1[1]+3))]),
                       max(data$DATA.12[c((data$PLOC.1[1]-3):(data$PLOC.1[1]+3))]))
                      )*58,digits=0),
        unlist(lapply(c(2:(length(data$PLOC.1)-1)),function(x){round((1-noise(max(data$DATA.9[c((data$PLOC.1[x]-3):(data$PLOC.1[x]+3))]),
                                                                              max(data$DATA.10[c((data$PLOC.1[x]-3):(data$PLOC.1[x]+3))]),
                                                                              max(data$DATA.11[c((data$PLOC.1[x]-3):(data$PLOC.1[x]+3))]),
                                                                              max(data$DATA.12[c((data$PLOC.1[x]-3):(data$PLOC.1[x]+3))]))
                                                                             )*58,digits=0)
                                                              }
                     )
              ),
        round((1-noise(max(data$DATA.9[c((data$PLOC.1[to]-3):(data$PLOC.1[to]+3))]),
                       max(data$DATA.10[c((data$PLOC.1[to]-3):(data$PLOC.1[to]+3))]),
                       max(data$DATA.11[c((data$PLOC.1[to]-3):(data$PLOC.1[to]+3))]),
                       max(data$DATA.12[c((data$PLOC.1[to]-3):(data$PLOC.1[to]+3))]))
                      )*58,digits=0)
    )
    return(qual)
}

align_intens_calls <- function(calls_fwd,intens_fwd,calls_rev=NULL,intens_rev=NULL){
    calls_intens_fwd <- numeric(length(calls_fwd))
    res_fwd <- lapply(1:4,function(x) {
        calls_intens_fwd[which(calls_fwd != "-")] <- intens_fwd[[x]]
        return(calls_intens_fwd)
        })
    if(!is.null(calls_rev)){
        calls_intens_rev <- numeric(length(calls_rev))
        res_rev <- lapply(4:1,function(x) {
            calls_intens_rev[which(calls_rev != "-")] <- rev(intens_rev[[x]])
            return(calls_intens_rev)
        })
        return(data.table(iA_fwd = res_fwd[[1]],iC_fwd = res_fwd[[2]],iG_fwd = res_fwd[[3]],iT_fwd = res_fwd[[4]],iA_rev = res_rev[[1]],iC_rev = res_rev[[2]],iG_rev = res_rev[[3]],iT_rev = res_rev[[4]]))
    }
    return(data.table(iA_fwd = res_fwd[[1]],iC_fwd = res_fwd[[2]],iG_fwd = res_fwd[[3]],iT_fwd = res_fwd[[4]]))
}

generate_ref <-function(user_seq,glassed_ref){
    #dummy return in case there is no reference 
    if(gsub("data/refs/","",gsub(".glassed.intrex.fasta","",glassed_ref))=="-"){
        return(list(data.table(exon_intron=NA,id=c(1:length(user_seq)),gen_coord=NA),data.table(id=numeric(0)),numeric(0)))
    }
    ref_aln_cover_init <- 0.75 # 0.85
    ref_aln_score_init <- 50 # 50
    ref_aln_lngth_gaps <- 50 # 50
    ref_aln_cover_gaps <- 0.75 # 0.6
    ref_aln_cover_bits <- 0.75 # 0.6
    ref_aln_score_bits <- 10 # 10

    cores         <- 2
    multiply_covered <- list()
#     user_seq <- gsub("[^ACGT]","N",user_seq)

    refs   <- readLines(glassed_ref)
    ref_info <- gsub(">ref_","",perl = T,refs[seq(1,length(refs),2)])
    ref_info <- strsplit(ref_info,split = "_")
    ref_names <- sapply(ref_info,function(x) x[1])
    ref_start <- as.numeric(sapply(ref_info,function(x) x[3]))
    ref_end <- as.numeric(sapply(ref_info,function(x) x[4]))

    align_split <- F
    
    if(align_split){
        refs <- DNAStringSet(refs[seq(2,length(refs),2)])
        align <- get_alignment(refs,user_seq,cores)
        OK_align <- which(align$score / nchar(refs) >= ref_aln_cover_init | align$score >= ref_aln_score_init) # OK_align <- which(align$score / nchar(refs) > 0.85 | align$score > 50)
        
        seq_coverage <- logical(nchar(user_seq))
        for(index in OK_align[order((align$score / nchar(refs))[OK_align],align$score[OK_align],decreasing = T)]){
            if(any(seq_coverage[(align$start[index] + 1):align$end[index]])) {
                indeces <- which(!seq_coverage[(align$start[index] + 1):align$end[index]]) + align$start[index]
                if(length(indeces) >= ref_aln_lngth_gaps | length(indeces) / (align$end[index] - align$start[index]) >= ref_aln_cover_gaps) { # if(length(indeces > 50) && length(indeces)/ (align$end[index] - align$start[index]) > 0.6 ) {
                    align$starts[index] <- min(indeces) - 1
                    align$ends[index] <- max(indeces)
                    ref_indeces <- which(!seq_coverage[(align$start[index] + 1):align$end[index]]) + align$ref_start[index] - 1
                    align$ref_start[index] <- min(ref_indeces)
                    align$ref_end[index] <- max(ref_indeces)
                    seq_coverage[indeces] <- T
                } else {
                    OK_align <- OK_align[-which(index == OK_align)]
                }
                multiply_covered[[ref_names[index]]] <- c(align$start[index],align$end[index],align$score[index])
            } else {
                seq_coverage[(align$start[index] + 1):align$end[index]] <- T
            }
        }
        
        all_aligns <- seq_along(align$score)
        for(index in all_aligns[order(align$score[all_aligns],decreasing = T)]){
            if(!any(seq_coverage[(align$start[index] + 1):align$end[index]]) && align$score[index] / (abs(align$end[index] - align$start[index])) >= ref_aln_cover_bits && align$score[index] >= ref_aln_score_bits) {
                seq_coverage[(align$start[index] + 1):align$end[index]] <- T
                OK_align <- c(OK_align,index)
            }
        }
        
        user_seq_vs_genome <- rbindlist(lapply(seq_along(OK_align),function(x) data.table(exon_intron = ref_names[OK_align][x]
                                                                                          ,id = sort(c((align$start[OK_align][x] + 1):align$end[OK_align][x],add_insert(align$diffs[type == "I" & id == OK_align[x]]) - 1))
                                                                                          ,gen_coord = get_coord(align$start[OK_align][x],align$ref_start[OK_align][x],align$ref_end[OK_align][x],ref_start[OK_align][x],ref_end[OK_align][x],align$diffs[type == "D" & id == OK_align[x]]))))
        
        #align$diffs[type == "D",t_pos := t_pos + 1]
        align$diffs <- align$diffs[id %in% OK_align,]
    } else {
        
        user_seq_vs_genome <- rbindlist(lapply(seq_along(ref_names),function(x) data.table(exon_intron = ref_names[x],gen_coord_abs = ref_start[x]:ref_end[x])))
        ref <- toupper(paste(refs[seq(2,length(refs),2)],collapse = ""))
        align <- get_alignment(ref,user_seq,cores,"local-global")

        gen_coord <- get_coord(align$start,align$ref_start,align$ref_end,max(user_seq_vs_genome$gen_coord),min(user_seq_vs_genome$gen_coord),align$diffs[type == "D"])
        #user_seq_vs_genome <- user_seq_vs_genome[align$ref_starts:align$ref_ends,]
        
        #merge does not keep the coorect order
        #user_seq_vs_genome <- merge(data.table(gen_coord = gen_coord,gen_coord_abs = ceiling(gen_coord)),user_seq_vs_genome,by = "gen_coord_abs")
        dt<- data.table(gen_coord = gen_coord,gen_coord_abs = ceiling(gen_coord))
        user_seq_vs_genome <- join(dt,user_seq_vs_genome,by="gen_coord_abs")
        user_seq_vs_genome[,gen_coord_abs := NULL]
        #user_seq_vs_genome[,id := sort(c((align$start + 1):align$end,add_insert(align$diffs[type == "I"]) - 1))]
        user_seq_vs_genome$id <-  sort(c((align$start + 1):align$end,add_insert(align$diffs[type == "I"]) - 1))
        setcolorder(user_seq_vs_genome,c(1,3,2))
    }
    
    return(list(user_seq_vs_genome,align$diffs,add_insert(align$diffs[type == "I"]) - 1))
}

add_insert <- function(diffs){
    if(nrow(diffs) > 0) return(unlist(lapply(1:nrow(diffs),function(x) diffs[["t_pos"]][x] + (1:nchar(diffs[["replace"]][x]))/100 )))
    else return(numeric())
}

get_coord <- function(seq_start,al_start,al_end,ref_start,ref_end,diffs){
    coord <- (ref_start:ref_end)[al_start:al_end]
    if(nrow(diffs) > 0){
        diffs[,t_pos := t_pos - seq_start]
        coord_out <- numeric(length(coord) + nrow(diffs))
        coord_out[-diffs[["t_pos"]]] <- coord
        for(index in which(coord_out == 0)) coord_out[index] <- coord_out[index - 1] - 0.01
        return(coord_out)
    } else return(coord)
}

get_fwd_rev_align <- function(fwd_seq,rev_seq,fwd_qual,rev_qual){
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

    #create quality scores corresponding with alignments
    fwd_split_qual <- rep(0,cons_length)
    rev_split_qual <- fwd_split_qual
    fwd_split_qual[which(fwd_split != "-")] <- fwd_qual
    rev_split_qual[which(rev_split != "-")] <- rev(rev_qual)

    fwd_offset <- start(pattern(align))
    rev_offset <- start(subject(align))

    #construct consensus /higher quality wins/
    cons_split <- fwd_split
    cons_split[which(fwd_split_qual < rev_split_qual)] <- rev_split[which(fwd_split_qual < rev_split_qual)]
    return(list(cons_split,fwd_split,rev_split,fwd_split_qual,rev_split_qual,fwd_offset,rev_offset))
}

get_alignment <- function(data,user_seq,cores,type = "local"){
    insert_weights <- c(1,4,8,12.6,16.3,25,30,rep(50,100))
    if(length(data) < 1) return(NULL)
    inserts_universal <- F
    if(length(data) <= cores * 2) cores <- 1
    splitid <- sort(seq_along(data) %% cores)

    a <- lapply(1:cores - 1L, function(x) {
        align <- pairwiseAlignment(pattern = data[which(splitid == x)], subject = user_seq,type = type,substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
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
get_intensities <- function(data,data_rev,calls,deletions=NULL,norm=FALSE,single_rev=FALSE) {
    #abi file documentation http://www.bioconductor.org/packages/release/bioc/vignettes/sangerseqR/inst/doc/sangerseq_walkthrough.pdf

    rev <- !(is.null(data_rev))

    intens <- data.table(data[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="A")),sep="")][[1]],
                         data[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="C")),sep="")][[1]],
                         data[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="G")),sep="")][[1]],
                         data[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="T")),sep="")][[1]])

    if(rev) intens_rev <- data.table(data_rev[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="A")),sep="")][[1]],
                                     data_rev[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="C")),sep="")][[1]],
                                     data_rev[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="G")),sep="")][[1]],
                                     data_rev[paste0("DATA.",toString(8+which(strsplit(data$FWO,split='')[[1]]=="T")),sep="")][[1]])
    else    intens_rev <- NULL

    #!!out of date
#     if(norm==TRUE){
#         #first and last 500 points
#         f_intens_start <- data.table()
#         f_intens_end   <- data.table()
#         last <- nrow(intens)
#         for(i in 1:499){
#             f_intens_start <- rbind(f_intens_start,as.list(apply(intens[1:(i+499)],2,function (x){sum(x)/(499+i)})))
#             f_intens_end   <- rbind(f_intens_end,as.list(apply(intens[(last-998+i):last],2,function(x){sum(x)/(998-i)})))
#         }
#         #1000 rolling mean
#         f_intens_mid <- data.table(rollmean(intens,k=999))
#         f_intens <- rbind(f_intens_start,f_intens_mid,f_intens_end)
#         intens<-intens/f_intens
#     }

    # clipping the end of forward chromatogram after last call
    if(single_rev){ #single reverse abi will be flipped over
        offset <- data$PLOC.1[1]
        #if(nrow(intens) > max(data$PLOC.1)){
        #    intens <- intens[1:max(data$PLOC.1)+offset]
        #}
    }else{
        if(nrow(intens) > max(data$PLOC.1)){
            intens <- intens[1:max(data$PLOC.1)]
        }
    }

    # clipping the end of reverse chromatogram after last call + length of first trace peak in fwd so we get the same offset once we flip it over
    if(rev){
        offset <- data$PLOC.1[1]
        if(nrow(intens_rev) > max(data_rev$PLOC.1)+offset){
            intens_rev <- intens_rev[1:max(data_rev$PLOC.1)+offset]
        }
    }
    peak_res <- 11
    intens <- normalize_peak_width(intens,data$PLOC.1,peak_res)
    calls  <- calls[,trace_peak:=rescale_call_positions(data$PLOC.1[1],nrow(calls),peak_res)]
    if(single_rev){
        intens <- setnames(data.table(intens),c("T","G","C","A"))
        deletions <- calls[call=="-"]
        intens <- intens[1:(max(calls$trace_peak) - (nrow(deletions))*peak_res + data$PLOC.1[1])]
        for (j in seq_len(ncol(intens)))
            set(intens,which(is.na(intens[[j]])),j,0)
        intens <- intens[nrow(intens):1,]
    } else {
        intens <- setnames(data.table(intens),c("A","C","G","T"))
    }
    intens[,id:=c(1:nrow(intens))]
    setkey(intens,id)
    if(rev){
        intens_rev <- normalize_peak_width(intens_rev,data_rev$PLOC.1,peak_res)
        intens_rev <- setnames(data.table(intens_rev),c("T","G","C","A"))
        intens_rev <- intens_rev[nrow(intens_rev):1]
        calls      <- calls[,trace_peak_rev:=rescale_call_positions(data_rev$PLOC.1[1],nrow(calls),peak_res)]
        intens_rev[,id:=c(1:nrow(intens_rev))]
        setkey(intens_rev,id)
    }

    deletions <- calls[call=="-"][,id]
    if(rev) deletions_rev <- calls[call_rev=="-"][,id]
    else deletions_rev <- list()
    if(length(deletions)!=0){

        del_pos <- 0
        del_pos <- calls[id %in% deletions][,trace_peak]
        rep <- 0
        for(i in c(1:length(del_pos))){
            pos <- del_pos[i]
            for(i in 1:12){intens<-rbind(intens,list(A=0,C=0,G=0,T=0,id=(pos-6 +rep+ i/100)))}
            rep <- rep - 12
        }
        setkey(intens,id)
        intens[,id:=c(1:nrow(intens))]
        setkey(intens,id)
    }
    if(length(deletions_rev)!=0){
        del_pos <- 0
        del_pos <- calls[id %in% deletions_rev][,trace_peak_rev]
        rep <- 0
        for(i in c(1:length(del_pos))){
            pos <- del_pos[i]
            for(i in 1:12){intens_rev<-rbind(intens_rev,list(A=0,C=0,G=0,T=0,id=(pos-6 +rep+ i/100)))}
            rep <- rep - 12
        }
        setkey(intens_rev,id)
        intens_rev[,id:=c(1:nrow(intens_rev))]
        setkey(intens_rev,id)
    }

    calls <- cbind(calls,get_intens(intens,calls[["trace_peak"]]))
    if(rev) calls <- cbind(calls,get_intens(intens_rev,calls[["trace_peak"]],"rev"))

    return(list(intens=intens,intens_rev=intens_rev,calls=calls))
}

get_intens <-function(intens,vec,tag = "fwd"){
    names <- c("A","C","G","T")
    #if(single_rev) names <- c("T","C","G","A")
    res <- lapply(names,function(x) get_single_intens(intens[[x]],vec))

    if(tag == "fwd") return(data.table(iA_fwd = res[[1]],iC_fwd = res[[2]],iG_fwd = res[[3]],iT_fwd = res[[4]]))
    else return(data.table(iA_rev = res[[1]],iC_rev = res[[2]],iG_rev = res[[3]],iT_rev = res[[4]]))
}

get_single_intens <- function(intens,call_vec){
    res <- sapply(call_vec,function(x) get_window_intens(intens[max((x-5),0):min((x+5),length(intens))]))
    res[is.na(res)] <- 0
    return(res)
}

get_window_intens <- function(intens){
    min_neg <- which((diff(intens) < 0))
    min_pos <- which((diff(intens) > 0))
    suppressWarnings(point <- which(min_neg > min(min_pos)))
    if(length(min_neg) == 0 || length(min_pos) == 0 || length(point) == 0 ) return(min(intens))
    return(intens[min_neg[min(point)]])
}

#if a mutation is around 50 % the variant peak may be sometimes higher than the reference peak
#then this function swaps the bases and percentages of the sample and mutation peak
adjust_ref_mut <- function(calls,intens_rev){
    vf <- calls[(mut_peak_base_fwd == reference ) & (abs(mut_peak_pct_fwd - sample_peak_pct_fwd) < 40 ), ]$id
    if (length(vf)>=1){
        calls[id %in% vf,`:=`(swp_base = mut_peak_base_fwd,swp_pct = mut_peak_pct_fwd)]
        calls[id %in% vf,`:=`(mut_peak_base_fwd = sample_peak_base_fwd,mut_peak_pct_fwd = sample_peak_pct_fwd)]
        calls[id %in% vf,`:=`(sample_peak_base_fwd =  swp_base,sample_peak_pct_fwd = swp_pct)]
        calls[id %in% vf,call := sample_peak_base_fwd]
        calls[,`:=`(swp_base = NULL,swp_pct = NULL)]
    }
    if(!is.null(intens_rev)){
        vr <- calls[(mut_peak_base_rev == reference ) & (abs(mut_peak_pct_rev - sample_peak_pct_rev) < 40 ), ]$id
        if (length(vr)>=1){
            calls[id %in% vr,`:=`(swp_base = mut_peak_base_rev,swp_pct = mut_peak_pct_rev)]
            calls[id %in% vr,`:=`(mut_peak_base_rev = sample_peak_base_rev,mut_peak_pct_rev = sample_peak_pct_rev)]
            calls[id %in% vr,`:=`(sample_peak_base_rev =  swp_base,sample_peak_pct_rev = swp_pct)]
            calls[id %in% vf,call_rev := sample_peak_base_rev]
            calls[,`:=`(swp_base = NULL,swp_pct = NULL)]
        }
    }
    return(calls)
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
hardtrim <- function(data,limit){
    last_int <- max(last(which(data$DATA.12 > limit)),last(which(data$DATA.11 > limit)),last(which(data$DATA.10 > limit)),last(which(data$DATA.9 > limit)))
    last_seq <- last(which(data$PLOC.1 < last_int))
    newdata = list()
    newdata$PLOC.1 <- data$PLOC.1[1:last_seq]
    newdata$PCON.1 <- data$PCON.1[1:last_seq]
    newdata$PCON.2 <- data$PCON.2[1:last_seq]
    newdata$PBAS.1 <- substr(data$PBAS.1,1,last_seq)
    newdata$DATA.9 <- data$DATA.9[1:last_int]
    newdata$DATA.10 <- data$DATA.10[1:last_int]
    newdata$DATA.11 <- data$DATA.11[1:last_int]
    newdata$DATA.12 <- data$DATA.12[1:last_int]
    newdata$FWO <- data$FWO
    
    return(newdata)

}
