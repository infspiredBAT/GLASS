library(zoo)
library(seqinr)
# the function extracts the signal intesities for each channel and returns it formatted for the javascript chromatograph
get_intensities <- function(data) {
  
    #abi file documentation http://www.bioconductor.org/packages/release/bioc/vignettes/sangerseqR/inst/doc/sangerseq_walkthrough.pdf
    
    abi_data <- data
    
    
    ins <- data.table(abi_data$DATA.9,abi_data$DATA.10,abi_data$DATA.11,abi_data$DATA.12)
    fwo <- abi_data$FWO
    setnames(ins,c(substring(fwo,1,1),substring(fwo,2,2),substring(fwo,3,3),substring(fwo,4,4)))
    return(ins)   
}

#function extracting 
get_call_data <- function(data){
  
    #TO DO if(length(data$PLOC.1)<=length(data$PBAS.1)){}
    qual      <- data$PCON.2
    rm7qual   <- rollmean(qual,k=7)
    rm7qualext<- data.frame(c(rep(rm7qual[1],3),rm7qual,rep(rm7qual[length(qual)-6],3)))  #cbind returns data frame if atleast one argument is a data frame
    res       <- generate_ref(data)
    seq_trim  <- str_split(data$PBAS.2,pattern="")[[1]]
    seq_trim[rm7qualext<12|qual<10]<-"low_qual"
    #save send to Vojta
    #+ fasta files with refs
    call      <- cbind(c(1:length(data$PLOC.2)),data$PLOC.2,
                       str_split(data$PBAS.2,pattern="")[[1]],
                       seq_trim,qual,rm7qualext,res$g_ref,stringsAsFactors = FALSE)
    colnames(call)<-c("id","trace_peak","call","seq_trim",
                      "quality","7rollmean_qual","reference",
                      "gen_coord","exon_intron")
    
    # changing the sequence coordinates to intensities coordinates for the brush tool
    res$helperdat$helper_intrex$start <- lapply(res$helperdat$helper_intrex$start,function (x) {call$trace_peak[x]})
    res$helperdat$helper_intrex$end   <- lapply(res$helperdat$helper_intrex$end  ,function (x) {call$trace_peak[x]})
    return(list(call=call,helperdat=res$helperdat))
    
}

generate_ref <-function(data){
    print(getwd())
    seq    <- data$PBAS.2
    refs   <- read.fasta("../../data/ref_ex_in.fa",as.string=TRUE)
    g_ref  <- c(rep("",nchar(data$PBAS.2)))     #generated reference
    coord  <- c(rep(NA,nchar(data$PBAS.2)))
    intrex <- c(rep(NA,nchar(data$PBAS.2)))
    helper_intrex   <- data.frame(attr=character(),start=numeric(),end=numeric()) #data for showing introns/exons in brush, maybe there is a cleverer way to do this
    #STEP 1.  
    for(ref in refs){
        m<-matchPattern(pattern = toupper(ref[1]),
                        subject = seq,
                        max.mismatch = 10,
                        min.mismatch=0,
                        with.indels = TRUE)
        
        if(length(m)!=0){
            if(m@ranges@width==
               as.integer(strsplit(attr(ref,"Annot"),split='_')[[1]][5])){ #if there no indels
                if(nmismatch(pattern=toupper(ref[1]),m)>0){                #if there are variations we must make sure
                    #print(mismatch(pattern=toupper(ref[1]),m))             #to assign ref
                    g_ref[start(m):end(m)]  <- strsplit(toupper(ref),"")[[1]]    
                }else{ #no variations                                                                   
                    g_ref[start(m):end(m)]  <- suppressWarnings(as.matrix(m))
                }
                coord[start(m):end(m)]  <- (strsplit(attr(ref,"name"),"_")[[1]][3]:
                                            strsplit(attr(ref,"name"),"_")[[1]][4])
                intrex[start(m):end(m)] <- strsplit(attr(ref,"name"),"_")[[1]][2]                  
                helper_intrex <- rbind(helper_intrex,
                                     data.frame(attr=strsplit(attr(ref,"name"),"_")[[1]][2],
                                     start=start(m),end=end(m)))
            }else{
              #deal with indels
            }
        }
    }
    max_y <- max(data$DATA.9,data$DATA.10,data$DATA.11,data$DATA.12)   #max_y prevoiusly calculeted in javascript, faster in R
    return(list(g_ref=cbind(g_ref,coord,intrex),helperdat=list(helper_intrex=helper_intrex,max_y=max_y)))
}

