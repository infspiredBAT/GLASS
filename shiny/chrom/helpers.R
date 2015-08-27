library(zoo)
library(seqinr)
# the function extracts the signal intesities for each channel and returns it formatted for the javascript chromatograph
get_intensities <- function(data) {
  
    #abi file documentation http://www.bioconductor.org/packages/release/bioc/vignettes/sangerseqR/inst/doc/sangerseq_walkthrough.pdf
    
    abi_data <- data
    
    #TO DO intensities data will be used for data analysis (not only visualisations)
    #create better R structure leave the conversion for later steps before visualisation
    
    
    #initialize long lists of x's and y's for naming the coordinates
#     y9 <-setNames(as.list(abi_data$DATA.9),
#                  rep(c("y"),times=length(abi_data$DATA.9)))
#     x9 <-setNames(as.list(c(1:length(abi_data$DATA.9))),
#                  rep(c("x"),times=length(abi_data$DATA.9)))
#     y10<-setNames(as.list(abi_data$DATA.10),
#                   rep(c("y"),times=length(abi_data$DATA.10)))
#     x10<-setNames(as.list(c(1:length(abi_data$DATA.10))),
#                   rep(c("x"),times=length(abi_data$DATA.10)))
#     y11<-setNames(as.list(abi_data$DATA.11),
#                   rep(c("y"),times=length(abi_data$DATA.11)))
#     x11<-setNames(as.list(c(1:length(abi_data$DATA.11))),
#                   rep(c("x"),times=length(abi_data$DATA.11)))
#     y12<-setNames(as.list(abi_data$DATA.12),
#                   rep(c("y"),times=length(abi_data$DATA.12)))
#     x12<-setNames(as.list(c(1:length(abi_data$DATA.12))),
#                   rep(c("x"),times=length(abi_data$DATA.12)))
#     
#     
#     
#     xylist9  <- list()
#     xylist10 <- list()
#     xylist11 <- list()
#     xylist12 <- list()
#     
#     for(i in 1:length(abi_data$DATA.9)){xylist9[[i]]<-append(x9[i],y9[i])}
#     for(i in 1:length(abi_data$DATA.10)){xylist10[[i]]<-append(x10[i],y10[i])}
#     for(i in 1:length(abi_data$DATA.11)){xylist11[[i]]<-append(x11[i],y11[i])}
#     for(i in 1:length(abi_data$DATA.12)){xylist12[[i]]<-append(x12[i],y12[i])}
#     
#     #extract the order of nucleotide order from abi file FWO_.1
#     list9 <-append(substring(abi_data$FWO,1,1),list(xylist9))
#     list10<-append(substring(abi_data$FWO,2,2),list(xylist10))
#     list11<-append(substring(abi_data$FWO,3,3),list(xylist11))
#     list12<-append(substring(abi_data$FWO,4,4),list(xylist12))
#     
#     #alternative without y
#     alist9 <-append(substring(abi_data$FWO,1,1),list(x9))
#     alist10<-append(substring(abi_data$FWO,2,2),list(x10))
#     alist11<-append(substring(abi_data$FWO,3,3),list(x11))
#     alist12<-append(substring(abi_data$FWO,4,4),list(x12))
#     
#     list9names <-setNames(list9,c("name","data"))
#     list10names<-setNames(list10,c("name","data"))
#     list11names<-setNames(list11,c("name","data"))
#     list12names<-setNames(list12,c("name","data"))
#     
#     alist9names <-setNames(alist9,c("name","data"))
#     alist10names<-setNames(alist10,c("name","data"))
#     alist11names<-setNames(alist11,c("name","data"))
#     alist12names<-setNames(alist12,c("name","data"))
#     
#     listALL<-list(list9names,list10names,list11names,list12names)
#     alistALL<-list(alist9names,alist10names,alist11names,alist12names)
    
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

