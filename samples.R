samplesLoad <- function(s_files,output,g_files,alignTo,g_custom_ref){
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
                            if(! alignTo[j] == "Custom"){
                                ref <- read.fasta(paste0("data/refs/",alignTo[j],".glassed.intrex.fasta"))
                            }else {
                                ref <- read.fasta(g_custom_ref)
                            }
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
                        loaded <- rbind(loaded,list(FWD_name=s_files[i,]$name,FWD_file=s_files[i,]$datapath,REV_name="-",REV_file="-",REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls="",brush_fwd_start= 0,brush_fwd_end=0,brush_rev_start=0,brush_rev_end = 0,coding = '',protein = '',VAF='',dbSNP = '', dbSNP_id = ''))
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
                                loaded <- rbind(loaded,list(FWD_name="-",FWD_file="-",REV_name=s_files[i,]$name,REV_file=s_files[i,]$datapath,REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls="",brush_fwd_start = 0,brush_fwd_end=0,brush_rev_start=0,brush_rev_end = 0,coding = '',protein = '',VAF = '', dbSNP = '', dbSNP_id = ''))
                            }
                        }else{
                            loaded <- rbind(loaded,list(FWD_name="-",FWD_file="-",REV_name=s_files[i,]$name,REV_file=s_files[i,]$datapath,REF=ref_name,mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,status="new",calls="",brush_fwd_start = 0,brush_fwd_end=0,brush_rev_start=0,brush_rev_end = 0,coding = '',protein = '',VAF = '',dbSNP = '', dbSNP_id = ''))
                        }
                    }

                }
                else not_loaded <- c(not_loaded,s_files[i,]$name)
            }
        }
        return(list(not_loaded=not_loaded,loaded=loaded))
    })
}

get_gbk_info <- function(session,file){
    #  ret=system("dev/GLASS/ext/gb2tab.py -a 1000 -b 1000 -f 'mRNA,CDS' Desktop/tp53.gb",intern=TRUE)
    #call <- paste0(c("python ext/gb2tab.py -f 'mRNA' ",file$datapath),collapse = "")
    #mRNA_call <- system(call,intern=TRUE)
    mRNA=NULL
    #for(i in c(1:length(mRNA_call))){
    #    tab <- unlist(strsplit(mRNA_call[i],"\t"))
    #    transcript_id <- gsub('transcript_id=\\"',"",str_match(tab[[4]],'transcript_id=\"[A-Z,a-z,0-9,_,-,.]+')[,1])
    #    spliced_product <- gsub('spliced_product=\\"',"",str_match(tab[[4]],'spliced_product=\"[A-Z,a-z]+')[,1])
    #    m <-list(transcript_id=transcript_id,spliced_product=spliced_product)
    #    name <- paste0("mRNA",i,collapse = "")
    #    mRNA <-  c(mRNA,list(mRNA[[name]]<-m))
    #}

    CDS_call <- system2("python",c("ext/gb2tab.py","-f","'CDS'",file$datapath),stdout=TRUE)
    #CDS_call  <- system(call,intern=TRUE)
    CDS=NULL
    if (length(CDS_call)==0){
        stop(paste0("Unable to retreive information from GenBank file: ", file$name))
    }
    n=1
    #system2 chops the output
    CDS_merge = character()
    CDS_merge[n] <- CDS_call[1]
    if(length(CDS_call > 1)){
        for(i in c(2:c(length(CDS_call)))){
            CDSt <- paste0(CDS_merge[n],CDS_call[i])
            if(length(unlist(strsplit(CDSt[1],"\t")))<5){
                CDS_merge[n] <- CDSt
            }else{
                n=n+1
                CDS_merge[n]<-CDS_call[i]
            }
        }
    }

    for(i in c(1:length(CDS_merge))){
        tab <- unlist(strsplit(CDS_merge[i],"\t"))
        transcript_id <- gsub('transcript_id=\\"',"",str_match(tab[[4]],'transcript_id=\"[A-Z,a-z,0-9,_,-,.]+')[,1])
        product <- gsub('product=\\"',"",str_match(tab[[4]],'product=\"[A-Z,a-z, ,0-9,_,-]+')[,1])
        translation <- gsub('translation=\\"',"",str_match(tab[[4]],'translation=\"[A-Z]+')[,1])
        gene <- gsub('gene=\\"',"",str_match(tab[[4]],'gene=\"[A-Z,a-z, ,0-9,_,-]+')[,1])
        m <-list(transcript_id=transcript_id,product=product,tabi=i,gene=gene,translation=translation)
        name <- paste0("CDS",i,collapse = "")
        CDS <-  c(CDS,list(CDS[[name]]<-m))
    }
    return(list("mRNA"=mRNA,"CDS"=CDS))
}

process_gbk <- function(session,file,ind){
    #INIT
    input_orient = "+"
    input_gene_start = 1
    input_chrom = "UN"
    input_gene_name = "UN"

    incProgress(1/10,message=NULL)

    CDS_call <- system2("python",c("ext/gb2tab.py","-f","'CDS'",file$datapath),stdout=TRUE)
    #CDS_call  <- system(call,intern=TRUE)
    CDS=NULL
    if (length(CDS_call)==0){
        stop(paste0("Unable to retreive information from GenBank file: ", file$name))
    }
    n=1
    #system2 chops the output
    CDS_merge = character()
    CDS_merge[n] <- CDS_call[1]
    if(length(CDS_call > 1)){
        for(i in c(2:c(length(CDS_call)))){
            CDSt <- paste0(CDS_merge[n],CDS_call[i])
            if(length(unlist(strsplit(CDSt[1],"\t")))<5){
                CDS_merge[n] <- CDSt
            }else{
                n=n+1
                CDS_merge[n]<-CDS_call[i]
            }
        }
    }

    #extract info from the genbank file

    tab1 <- unlist(strsplit(CDS_merge[ind],"\t"))

    incProgress(1/10,message=NULL)

    #coordinates <- str_match(tab1[[4]],'/GenBank.* REGION: [0-9]+..[0-9]+')[,1]
    #input_chrom <- gsub("NC_0+","",str_match(coordinates,"NC_0+[1-9,X]+"))      #might not work for chr X,Y and MT
    #input_orient <- gsub("strand=\\\"","",str_match(tab1[[4]],"strand=\\\"."))
    #if(input_orient=="+"){
    #    input_gene_start <- as.numeric(gsub("REGION: ","",str_match(coordinates,"REGION: [0-9]+")))
    #}else{
    #    input_gene_start <- as.numeric(unlist(str_match_all(coordinates,"[0-9]+"))[length(unlist(str_match_all(coordinates,"[0-9]+")))])
    #}
    #if(is.na(coordinates)){
    #    coordinates <- str_match(tab1[[4]],'/GenBank.* REGION: complement.[0-9]+..[0-9]+')[,1]
    #    input_orient <- "-"
    #    input_chrom <- gsub("NC_0+","",str_match(coordinates,"NC_0+[1-9]+"))
    #    input_gene_start <- as.numeric(unlist(str_match_all(coordinates,"[0-9]+"))[length(unlist(str_match_all(coordinates,"[0-9]+")))])
    #}

    exon <- data.table(which(strsplit(as.character(tab1[[3]]), '')[[1]]=='('),which(strsplit(as.character(tab1[[3]]), '')[[1]]==')'))
    setnames(exon,c("start","end"))
    exon[,id := seq_along(exon$start)]
    exon[,ex:=rep("exon",nrow(exon))]
    exon[,name := paste0(ex,id)]


    intron<-data.table(which(strsplit(as.character(tab1[[3]]), '')[[1]]=='D'),which(strsplit(as.character(tab1[[3]]), '')[[1]]=='A'))
    setnames(intron,c("start","end"))
    intron[,id := seq_along(intron$start)]
    intron[,ex:=rep("intron",nrow(intron))]
    intron[ ,`:=`( name =paste0(ex,id))]


    genedt<-rbind(exon,intron)
    setkey(genedt,start)
    genedt[,seq := substring(as.character(tab1[[2]]),start,end),by=1:nrow(genedt)]

    #strand specific

    if(input_orient == "+"){
        genedt[,start_chr := input_gene_start + start -1]
        genedt[,end_chr := input_gene_start + end -1]
    }else{
        genedt[,start_chr := input_gene_start - start +1]
        genedt[,end_chr := input_gene_start - end +1]
    }

    genedt[,len := end-start +1]
    genedt[,chr := rep(input_chrom,nrow(genedt))]
    custom_fasta <- tempfile()
    for(i in nrow(genedt)){
        writeLines(paste0(">ref_",genedt$name,"_",genedt$chr,"_",genedt$start_chr,"_",genedt$end_chr,"_",genedt$end-genedt$start,"\n",genedt$seq),custom_fasta)
    }
    #close(con) #fasta
    incProgress(1/10,message=NULL)

    cod_table <- rbindlist(
        lapply(1:nrow(genedt),
               function(x) data.table(codon = genedt$ex[x],
                                      gen_coord = genedt$start_chr[x]:genedt$end_chr[x],
                                      coding_seq = as.character(0),
                                      intrex =  as.character(genedt$id[x]),
                                      ord_in_cod = as.integer(0),AA = "",
                                      seq = unlist(strsplit(genedt$seq[x],split = "")))
               )
        )

    #ret=system("dev/GLASS/ext/gb2tab.py -a 20000 -b 20000 -f 'mRNA,CDS' Desktop/tp53.gb",intern=TRUE)
    #tab2  <- unlist(strsplit(ret[3],"\t"))
    #vec   <- strsplit(tab2[2],"")
    #vecb  <- lapply(vec,function(x){x==toupper(x)})[[1]]

    tab2 <- tab1
    trans <- gsub('/translation=\\"',"",str_match(tab2[[4]],'/translation=\"[A-Z]+')[,1])

    coding<-data.table(AA =c(unlist(strsplit(trans,split="")),"*"),codon = 1:length(c(unlist(strsplit(trans,split="")),"*")))
    coding <- rbindlist(lapply(1:nrow(coding),function(x) data.table(AA = rep(coding$AA[x],3),codon=as.character(rep(coding$codon[x],3)),ord_in_cod= as.integer(c(1,2,3)))))

    #remove non coding exons
    # cod_table <- cbind(cod_table, vecb)

    # cod_table[!vecb][codon=="exon"]$codon = "non-coding_exon_seq"
    incProgress(1/10,message=NULL)

    cod_table[codon=="exon"]$AA = coding$AA
    cod_table[codon=="exon"]$ord_in_cod = coding$ord_in_cod
    cod_table[codon=="exon",]$coding_seq = as.character(1:nrow(cod_table[codon == "exon",]))
    cod_table[codon=="exon"]$codon = coding$codon


    extra <- cod_table
    plus = extra[,max(as.numeric(coding_seq)),by=intrex][,V1]
    plus <- plus[1:length(plus)-1]
    minus <- plus + 1


    l1 <-unlist(lapply(1:length(plus),function(x) rep(plus[x],genedt[ex=="intron"]$len[x])))
    l2 <-unlist(lapply(1:length(minus),function(x) rep(minus[x],genedt[ex=="intron"]$len[x])))

    intr_cod <- data.table(plus = unlist(lapply(1:length(plus),function(x) rep(plus[x],genedt[ex=="intron"]$len[x]))),minus=unlist(lapply(1:length(minus),function(x) rep(minus[x],genedt[ex=="intron"]$len[x]))),p = unlist(lapply(1:length(plus), function(x) 1:genedt[ex=="intron"]$len[x] )),m = unlist(lapply(1:length(minus), function(x) genedt[ex=="intron"]$len[x]:1 )))

    #intron_name <- unlist(lapply(1:nrow(intr_cod), 
    #                                 function(x){
    #                                    if(intr_cod[x,]$p<intr_cod[x,]$m){ 
    #                                        return(paste0(intr_cod[x,]$plus,"+",intr_cod[x,]$p))
    #                                    } else {
    #                                        return(paste0(intr_cod[x,]$minus,"-",intr_cod[x,]$m))
    #                                    }
    #                                 }
    #                             )
    #                      )
    intr_cod[p<m,name:= paste0(plus,"+",p)]
    intr_cod[m<p,name:= paste0(minus,"-",m)]
    cod_table[codon=="intron"]$coding_seq = intr_cod$name

    incProgress(1/10,message=NULL)

    cod_table<- cod_table
    file_cod=paste0(input_gene_name,".glassed.codons.rdata")
    custom_cod <- tempfile()
    save(cod_table,file=custom_cod)

    incProgress(amount=1,message="Done")

    return(list(custom_cod=custom_cod,custom_fasta=custom_fasta))

}
