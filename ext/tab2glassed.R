#input = chromosome number, gene name, chromosome coordinate of first base, orientation
#edit manually
input_orient = "plus"  
input_gene_start = 13049414
input_chrom = 19
input_gene_name = "calreticulin"


tab <- read.table("CALR_DNA.tab",sep="\t")
genedt<-data.table(which(strsplit(as.character(tab[[3]]), '')[[1]]=='('),which(strsplit(as.character(tab[[3]]), '')[[1]]==')'))
setnames(genedt,c("start","end"))
genedt[,id := seq_along(genedt$start)]
genedt[,ex:=rep("exon",nrow(genedt))]
genedt[,name := paste0(ex,id)]


intron<-data.table(which(strsplit(as.character(tab[[3]]), '')[[1]]=='D'),which(strsplit(as.character(tab[[3]]), '')[[1]]=='A'))
setnames(intron,c("start","end"))
intron[,id := seq_along(intron$start)]
intron[,ex:=rep("intron",nrow(intron))]
intron[ ,`:=`( name =paste0(ex,id))]


genedt<-rbind(genedt,intron)
setkey(genedt,start)
genedt[,seq := substring(as.character(tab[[2]]),start,end),by=1:nrow(genedt)]

#strand specific

if(input_orient == "plus"){
    genedt[,start_chr := input_gene_start + start -1]
    genedt[,end_chr := input_gene_start + end -1]
}else{
    genedt[,start_chr := input_gene_start - start +1]
    genedt[,end_chr := input_gene_start - end +1]
}
genedt[,len := end-start +1]
genedt[,chr := rep(input_chrom,nrow(genedt))]
con<-file(paste0("CALR",".glassed.intrex.fasta"))
for(i in nrow(genedt)){writeLines(paste0(">ref_",genedt$name,"_",genedt$chr,"_",genedt$start_chr,"_",genedt$end_chr,"_",genedt$end-genedt$start,"\n",genedt$seq),con)}
close(con) #fasta



cod_table <- rbindlist(lapply(1:nrow(genedt),function(x) data.table(codon = genedt$ex[x],gen_coord = genedt$start_chr[x]:genedt$end_chr[x],coding_seq = as.character(0),intrex =  as.character(genedt$id[x]),ord_in_cod = as.integer(0),AA = "",seq = unlist(strsplit(genedt$seq[x],split = "")))))

#numbers differ in different genebank files grep for smth??
trans <- unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(tab[[4]]),";"))[5]," "))[11],"="))[9]
coding<-data.table(AA =c(unlist(strsplit(trans,split="")),"*"),codon = 1:length(c(unlist(strsplit(trans,split="")),"*")))
coding <- rbindlist(lapply(1:nrow(coding),function(x) data.table(AA = rep(coding$AA[x],3),codon=as.character(rep(coding$codon[x],3)),ord_in_cod= as.integer(c(1,2,3)))))

cod_table[codon=="exon"]$AA = coding$AA
cod_table[codon=="exon"]$ord_in_cod = coding$ord_in_cod
cod_table[codon=="exon",]$coding_seq = as.character(1:nrow(cod_table[codon == "exon",]))
cod_table[codon=="exon"]$codon = coding$codon


extra <- cod_table
plus = extra[,max(as.numeric(coding_seq)),by=intrex][,V1]
plus <- plus[1:length(plus)-1]
minus <- plus + 1


l1 <-unlist(lapply(1:length(plus),function(x) rep(plus[x],genedt[ex=="intron"]$len[x])))
l2<-unlist(lapply(1:length(minus),function(x) rep(minus[x],genedt[ex=="intron"]$len[x])))

intr_cod <- data.table(plus = unlist(lapply(1:length(plus),function(x) rep(plus[x],genedt[ex=="intron"]$len[x]))),minus=unlist(lapply(1:length(minus),function(x) rep(minus[x],genedt[ex=="intron"]$len[x]))),p = unlist(lapply(1:length(plus), function(x) 1:genedt[ex=="intron"]$len[x] )),m = unlist(lapply(1:length(minus), function(x) genedt[ex=="intron"]$len[x]:1 )))

name <- unlist(lapply(1:nrow(intr_cod), function(x) {if(intr_cod[x,]$p<intr_cod[x,]$m){ return(paste0(intr_cod[x,]$plus,"+",intr_cod[x,]$p))} else {return(paste0(intr_cod[x,]$minus,"-",intr_cod[x,]$m))}}))
cod_table[codon=="intron"]$coding_seq = name

cod_table<- cod_table
file_cod=paste0(input_gene_name,".glassed.codons.rdata") 
save(cod_table,file=file_cod)
