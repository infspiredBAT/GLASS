samples_load <- function(s_files,output,g_files){
    not_loaded <- ""
    loaded <- data.table()
    withProgress(message="processing...",value=0,{
        sm <- matrix(c(1 ,-1 ,-1 ,-1 ,-1 ,0.5 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,0.1 ,0 ,-1 ,1 ,-1 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.5 ,-1 ,0.1 ,-1 ,0.1 ,0.1 ,0 ,-1 ,-1 ,1 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,-1 ,0.1 ,0 ,-1 ,-1 ,-1 ,1 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.1 ,0.1 ,0.1 ,-1 ,0 ,-1 ,-1 ,0.5 ,0.5 ,0.1 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.5 ,0.5 ,-1 ,-1 ,-1 ,0.1 ,0 ,0 ,0 ,0 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.1 ,0.5 ,-1 ,0.5 ,-1 ,0 ,0 ,0.1 ,-1 ,0 ,0 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-1 ,0.5 ,-1 ,0.5 ,0 ,0 ,-1 ,0.1 ,0 ,0 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.5 ,0.5 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,-1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.5 ,-1 ,-1 ,0.5 ,0 ,0 ,0 ,0 ,-1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0 ,0.1 ,0 ,0 ,0.1 ,0.1 ,0.1 ,-1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0.1 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,-1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1),15,15,dimnames = list(c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N")))

        for(i in 1:nrow(s_files)){
            if(s_files[i,]$name %in% g_files$FWD_name | s_files[i,]$name %in% g_files$REV_name){
                not_loaded <- c(not_loaded,s_files[i,]$name)
            }else{
                incProgress(0,detail = paste("file ",i," of ",nrow(s_files)))
                abi<-NULL
                tryCatch(
                    abi <- sangerseqR::read.abif(s_files[i,]$datapath)@data
                )
                if(!is.null(abi)) {
                    seq <- DNAString(gsub("\\*","N",abi$PBAS.1))
                    mas <- 30                   #Minimum alignment score (Guess)
                    score_bst <- 0
                    ref_name <- "unk"
                    rev <- FALSE
                    output$samples_load_info <- renderPrint({paste0("processing file",i,sep="")})
                    for(j in 1:length(g_refs_avail)){
                        ref <- read.fasta(paste0("data/refs/",g_refs_avail[j],".glassed.intrex.fasta"))
                        ref <- toupper(paste0(unlist(ref),collapse = ""))
                        pa <- pairwiseAlignment(pattern = ref, subject = seq,type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
                        score_fwd <- pa@score
                        pa <- pairwiseAlignment(pattern = ref, subject = reverseComplement(seq),type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
                        score_rev <- pa@score
                        if(score_fwd > score_bst| score_rev > score_bst){
                            if(score_fwd > mas & score_fwd > score_rev){
                                ref_name <- g_refs_avail[j]
                                rev <- FALSE
                            }else if(score_rev > mas & score_rev > score_bst){
                                ref_name <- g_refs_avail[j]
                                rev <- TRUE
                            }
                        }
                    }
                    incProgress(1/nrow(s_files))
                    #test fwd/rev/ which reference
                    if(!rev)
                        loaded <- rbind(loaded,list(FWD_name=s_files[i,]$name,FWD_file=s_files[i,]$datapath,REV_name="-",REV_file="-",REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls=""))
                    else
                        loaded <- rbind(loaded,list(FWD_name="-",FWD_file="-",REV_name=s_files[i,]$name,REV_file=s_files[i,]$datapath,REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls=""))
                }
                else not_loaded <- c(not_loaded,s_files[i,]$name)
            }
        }
        return(list(not_loaded=not_loaded,loaded=loaded))
    })
}

