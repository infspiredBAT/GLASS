samples_load <- function(s_files,output,g_files,alignTo){
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
                seq<-NULL
                tryCatch(
                    {abi <- sangerseqR::read.abif(s_files[i,]$datapath)@data
                    seq <- DNAString(gsub("[\\*,!]","N",abi$PBAS.1))},
                    error = function(e){output$files <- renderPrint(paste0("<pre>error while loading abi file : ",e$message,"</pre>" ))}
                )
                if(!is.null(abi)&&!is.null(seq)) {
                    #try catch
                    #seq <- DNAString(gsub("\\*","N",abi$PBAS.1))
                    mas <- 60                   #Minimum alignment score (guess)
                    score_bst <- 0
                    ref_name <- "-"
                    rev <- FALSE
                    output$samples_load_info <- renderPrint({paste0("processing file",i,sep="")})
                    if(is.null(alignTo)){
                        ref <- "-"
                    }else{
                        #find best matching reference
                        for(j in 1:length(alignTo)){
                            ref <- read.fasta(paste0("data/refs/",alignTo[j],".glassed.intrex.fasta"))
                            ref <- toupper(paste0(unlist(ref),collapse = ""))
                            pa <- pairwiseAlignment(pattern = ref, subject = seq,type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
                            score_fwd <- pa@score
                            #print(paste0("ref: ", alignTo[j]))
                            #print(paste0("bst: ",score_bst))
                            #print(paste0("fwd: ",pa@score))
                            pa <- pairwiseAlignment(pattern = ref, subject = reverseComplement(seq),type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
                            score_rev <- pa@score
                            #print(paste0("rev: ",pa@score))
                            
                            if(score_fwd > score_bst || score_rev > score_bst){
                                
                                if(score_fwd > score_rev){
                                    score_bst <- score_fwd
                                }else{
                                    score_bst <- score_rev
                                }
                                if(score_fwd > mas & score_fwd > score_rev){
                                    ref_name <- alignTo[j]
                                    rev <- FALSE
                                }else if(score_rev > mas & score_rev > score_fwd){
                                    ref_name <- alignTo[j]
                                    rev <- TRUE
                                }
                            }
                        }
                    }
                    incProgress(1/nrow(s_files))
                    #test fwd/rev/ which reference
                    if(!rev)
                        loaded <- rbind(loaded,list(FWD_name=s_files[i,]$name,FWD_file=s_files[i,]$datapath,REV_name="-",REV_file="-",REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls="",brush_fwd = 0,brush_rev = 0,coding = '',protein = '',VAF='',dbSNP = '', dbSNP_id = ''))
                    else{
                        #pair fwd and rev samples with matching names
                        #if exactly one fwd file matches the stripped name of given rev file and this fwd file has no pair then:
                        
                        
                        rev_name <- gsub("_[12][09][0-9][0-9]-[0-1][0-9]-[0123][0-9]_[0-1][0-9]-[0-5][0-9]-[0-5][0-9]","",s_files[i,]$name) #remove date
                        rev_name <- gsub(".abi",".ab1",rev_name) #normalize extension
                        rev_name <- gsub(".ab1","",rev_name)     #remove extension
                        if(substr(rev_name,nchar(rev_name),nchar(rev_name))=="R" || substr(rev_name,nchar(rev_name),nchar(rev_name))=="F")
                            rev_name <- substr(rev_name,1,nchar(rev_name)-1) #remove last letter "R"/"F" naming style
                        matches <- grep(rev_name,loaded$FWD_name)
                        if(length(matches)==1){
                            if(loaded[matches]$REV_name == "-"&& loaded[matches]$REF==ref_name)
                                loaded[matches,`:=`(REV_name=s_files[i,]$name,REV_file=s_files[i,]$datapath)]
                            else{
                                loaded <- rbind(loaded,list(FWD_name="-",FWD_file="-",REV_name=s_files[i,]$name,REV_file=s_files[i,]$datapath,REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls="",brush_fwd = 0,brush_rev = 0,coding = '',protein = '',VAF = '', dbSNP = '', dbSNP_id = ''))
                            }
                        }else{
                            loaded <- rbind(loaded,list(FWD_name="-",FWD_file="-",REV_name=s_files[i,]$name,REV_file=s_files[i,]$datapath,REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls="",brush_fwd = 0,brush_rev = 0,coding = '',protein = '',VAF = '',dbSNP = '', dbSNP_id = ''))
                        }
                    }
                        
                }
                else not_loaded <- c(not_loaded,s_files[i,]$name)
            }
        }
        
        
        return(list(not_loaded=not_loaded,loaded=loaded))
    })
}

