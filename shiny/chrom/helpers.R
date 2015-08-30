library(zoo)
library(seqinr)

# the function extracts the signal intesities for each channel and returns it formatted for the javascript chromatograph
get_intensities <- function(data,norm=FALSE) {
    #abi file documentation http://www.bioconductor.org/packages/release/bioc/vignettes/sangerseqR/inst/doc/sangerseq_walkthrough.pdf
    intens <- data.table(data$DATA.9,data$DATA.10,data$DATA.11,data$DATA.12)
    fwo <- data$FWO
    setnames(intens,c(substring(fwo,1,1),substring(fwo,2,2),substring(fwo,3,3),substring(fwo,4,4)))

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
    return(intens)
}

get_call_data <- function(data){
    #TO DO if(length(data$PLOC.1)<=length(data$PBAS.1)){}
    qual      <- data$PCON.2
    rm7qual   <- rollmean(qual,k=7)
    #rm7qualext<- data.frame(c(rep(rm7qual[1],3),rm7qual,rep(rm7qual[length(qual)-6],3)))  #cbind returns data frame if atleast one argument is a data frame
    res       <- generate_ref(data)
    call <- data.table(id = seq_along(data$PLOC.2)
                   ,trace_peak = data$PLOC.2
                   ,call = str_split(data$PBAS.2,pattern="")[[1]]
                   ,seq_trim = str_split(data$PBAS.2,pattern="")[[1]]
                   ,quality = qual
                   ,reference = str_split(data$PBAS.2,pattern="")[[1]])
    call[,rm7qual := c(quality[1:3],rollmean(quality,k=7),quality[(length(quality) - 2):length(quality)])]
    setkey(res[[2]],id)
    if(length(res[[3]]) > 0) {
        setkey(call,id)
        add <- call[as.integer(res[[3]]),]
        data.table::set(add,NULL,"reference",unlist(strsplit(res[[2]][type == "I"][["replace"]],"")))
        call <- rbind(call,add[,id := res[[3]]][,call := "-"][,seq_trim := "-"])
    }
    call <- merge(call,res[[1]],all.x = T,by = "id")
    data.table::set(call,which(call[["id"]] %in% res[[2]][type != "I"][["t_pos"]]),"reference",res[[2]][type != "I"][["replace"]])
    data.table::set(call,which(is.na(call[["gen_coord"]])),"reference","NA")
    data.table::set(call,which(call[["rm7qual"]] < 12 | call[["quality"]] < 10),"seq_trim","low_qual")

    # changing the sequence coordinates to intensities coordinates for the brush tool
    helperdat <- list()
    helperdat$helper_intrex <- list()
    helperdat$helper_intrex <- setnames(call[!is.na(exon_intron),list(min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","start","end"))
    return(list(call=call,helperdat=helperdat))
}

generate_ref <-function(data){
#    print(getwd())
    cores <- 2

    user_seq <- gsub("[^ACGT]","N",data$PBAS.2)
#    user_seq <- 'TGNAGGAGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCANATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGGTCAGGAGCCACTTGCCACCCTGCACACTGGCCTGCTGTGCCCCAGCCTCTGCTTGCCTCTGACCCCTGGGCCCACCTCTTACCGATTTCTTCCATACTACTACCCATCCACCTCTCATCACATCCCCGGCGGGGAATCTCCTTACTGCTCCCACTCAGTTTTCTTTTCTCTGGCTTTGGGACCTCTTAACCTGTGGCTTCTCCTCCACCTACCTGGAGCTGGAGCTTAGGCTCCAGAAAGGACAAGGGTGGTTGGGAGTAGATGGAGCCTGGTTTTTTAAATGGGACAGGTAGGACCTGATTTCCTTACTGCCTCTTGCTTCTCTTTTCCTATCCTGGTAGTGGTAATCTTCTGGGACGGAACAGCTTTGAGGTTCATGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGTAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGCCTTGGACTCAGGATGCCCAGGCTGGAAGGAGCCAGGGGGAGCAGGCTCACTCCAGCCACCTGAAGTCCAAAANAGGGTA'
    refs   <- readLines("../../data/ref_ex_in.fa")
    ref_info <- gsub(">ref_","",perl = T,refs[seq(1,length(refs),2)])
    ref_info <- strsplit(ref_info,split = "_")
    ref_names <- sapply(ref_info,function(x) x[1])
    ref_start <- as.numeric(sapply(ref_info,function(x) x[2]))
    ref_end <- as.numeric(sapply(ref_info,function(x) x[3]))
    refs <- DNAStringSet(refs[seq(2,length(refs),2)])
    align <- get_alignment(refs,user_seq,cores)
    OK_align <- which(align$score / nchar(refs) > 0.2)

    user_seq_vs_genome <- rbindlist(lapply(seq_along(OK_align),function(x) data.table(exon_intron = ref_names[OK_align][x]
		,id = sort(c((align$start[OK_align][x] + 1):align$end[OK_align][x],add_insert(align$diffs[type == "I" & id == OK_align[x]])))
		,gen_coord = get_coord(align$start[OK_align][x],align$ref_start[OK_align][x],align$ref_end[OK_align][x],ref_start[OK_align][x],ref_end[OK_align][x],align$diffs[type == "D" & id == OK_align[x]]))))

    align$diffs[type == "D",t_pos := t_pos + 1]
    return(list(user_seq_vs_genome,align$diffs[id %in% OK_align,],add_insert(align$diffs[type == "I"])))

# #     g_ref  <- c(rep("",nchar(data$PBAS.2)))     #generated reference
# #     coord  <- c(rep(NA,nchar(data$PBAS.2)))
# #     intrex <- c(rep(NA,nchar(data$PBAS.2)))
# #     helper_intrex   <- data.frame(attr=character(),start=numeric(),end=numeric()) #data for showing introns/exons in brush, maybe there is a cleverer way to do this
#     #STEP 1.
#     for(ref in refs){
#         m<-matchPattern(pattern = toupper(ref[1]),
#                         subject = user_seq,
#                         max.mismatch = 10,
#                         min.mismatch=0,
#                         with.indels = TRUE)
#
#         if(length(m)!=0){
#             if(m@ranges@width==
#                as.integer(strsplit(attr(ref,"Annot"),split='_')[[1]][5])){ #if there no indels
#                 if(nmismatch(pattern=toupper(ref[1]),m)>0){                #if there are variations we must make sure
#                     #print(mismatch(pattern=toupper(ref[1]),m))             #to assign ref
#                     g_ref[start(m):end(m)]  <- strsplit(toupper(ref),"")[[1]]
#                 }else{ #no variations
#                     g_ref[start(m):end(m)]  <- suppressWarnings(as.matrix(m))
#                 }
#                 coord[start(m):end(m)]  <- (strsplit(attr(ref,"name"),"_")[[1]][3]:
#                                             strsplit(attr(ref,"name"),"_")[[1]][4])
#                 intrex[start(m):end(m)] <- strsplit(attr(ref,"name"),"_")[[1]][2]
#                 helper_intrex <- rbind(helper_intrex,
#                                      data.frame(attr=strsplit(attr(ref,"name"),"_")[[1]][2],
#                                      start=start(m),end=end(m)))
#             }else{
#               #deal with indels
#             }
#         }
#     }
#     max_y <- max(data$DATA.9,data$DATA.10,data$DATA.11,data$DATA.12)   #max_y prevoiusly calculeted in javascript, faster in R
#     return(list(g_ref=cbind(g_ref,coord,intrex),helperdat=list(helper_intrex=helper_intrex,max_y=max_y)))
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

get_alignment <- function(data,user_seq,cores,type = "overlap"){
    insert_weights <- c(1,4,8,12.6,16.3,25,30,rep(50,100))
    if(length(data) < 1) return(NULL)
    inserts_universal <- F
    if(length(data) <= cores * 2) cores <- 1
    splitid <- sort(seq_along(data) %% cores)
    sm <- matrix(-1,5,5,dimnames = list(c("A","C","G","T","N"),c("A","C","G","T","N")))
    diag(sm) <- 1
    sm[,"N"] <- 1
    sm["N",] <- 1
    a <- lapply(1:cores - 1L, function(x) {
        align <- pairwiseAlignment(pattern = data[which(splitid == x)], subject = user_seq,type = "overlap",substitutionMatrix = sm,gapOpening = -6, gapExtension = -0.3)
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
