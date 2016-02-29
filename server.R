library(shiny)
library(sangerseqR)
library(xlsx)
source("procAbi.R")
source("helpers.R")
source("samples.R")


shinyServer(function(input,output,session) {
    #options(shiny.reactlog=TRUE)
    g_calls                 <- NULL             #annotated basecall data
    #makeReactiveBinding("g_calls")
    g_intens                <- NULL             #intensities file
    g_intens_rev            <- NULL             #optional reverse intensities file
    g_intrexdat             <- NULL             #intrex data used in graphs
    g_glassed_ref           <- NULL
    g_glassed_cod           <- NULL
    g_choices               <- NULL
    g_noisy_neighbors       <- NULL
    g_view                  <- NULL
    g_selected              <- NULL
    g_selected_goto_index   <- 0
    g_max_y                 <- NULL
    g_hetero_calls          <- 0
    g_hetero_indel_pid      <- 0
    g_hetero_ins_tab        <- NULL
    g_hetero_del_tab        <- NULL
    g_expected_het_indel    <- NULL
    g_minor_het_insertions  <- NULL
    load_id                 <- NULL
    g_stored_het_indels     <- list()
    files_info              <- ""
    g_indels_present        <- FALSE
    g_qual_present          <- FALSE
    g_not_loaded            <- ""
    g_refs_avail            <<- c("TP53","NOTCH1","ATM")
    g_files                 <- data.table(FWD_name=c("LowFreq_frameShiftF.ab1"),
                                          FWD_file=c("data/abis/eric/3low_freq_fsF.ab1"),
                                          REV_name=c("lowFreq_frameShiftR.ab1"),
                                          REV_file=c("data/abis/eric/3low_freq_fsR.ab1"),
                                          REF=c("TP53")
    )
#     get_file <- reactive({
#        # if (!is.null(input$select_file)) return(input$select_file$datapath)
#         if (!is.null(input$select_file)) return(input$select_file)
#         else return(NULL)
#     })
    ex <- NULL
    btn_counter <- 0
    
    #SAMPLE BROWSER STUFF in progress...
    
    output$samples_table <- DT::renderDataTable({
            if(!is.null(input$browser_files)){
                g_not_loaded <- ""
                loaded <- ""
                ret <- samples_load(input$browser_files,output)
                g_files <<- rbind(g_files,ret$loaded)
            }
            disabled <- rep(FALSE,nrow(g_files))
            #dsbl[1] <- "true"
            add_load_buttons    <- shinyInput(actionButton, 1:nrow(g_files), 'loadSample_', label = NULL, onclick = 'Shiny.onInputChange(\"goLoadSamples\",  this.id)',ico=rep("play",nrow(g_files)) )
            add_delete_buttons    <- shinyInput(actionButton, 1:nrow(g_files), 'delSample_', label = NULL, onclick = 'Shiny.onInputChange(\"goDeleteSamples\",  this.id)',ico=rep("close",nrow(g_files)),dsbl=disabled )
            add_reference_dropdown <- shinyInput(selectInput, 1:nrow(g_files), 'selectInput_',choices=c("TP53","NOTCH1","ATM"),label=NULL,width="100px")
            out<-cbind(g_files[,list("forward"=FWD_name,"reverse"=REV_name)],"reference"=add_reference_dropdown,delete=add_delete_buttons,load=add_load_buttons)
            #DT::datatable(out,selection = "none")
        },
        escape=FALSE,
        options=list("paging"=FALSE,"searching"=FALSE,"autoWidth"=FALSE),selection="none"
    )
    
    
#    change_reference <- observe ({
#        input$gene_of_interest
#        isolate({
#            g_glassed_ref <<- paste("data/refs/",input$gene_of_interest,".glassed.intrex.fasta",sep="")
#            g_glassed_cod <<- paste("data/refs/",input$gene_of_interest,".glassed.codons.rdata",sep="")
#            #g_files <<- paste0("Reference changed to ", input$gene_of_interest)
#            #if(class(loading_processed_files())[1] != "my_UI_exception"){
#                #output$files      <-  renderPrint({cat(g_files)})
#            #}
#        })
#        #return(list(glassed_ref=glassed_ref,glassed_cod=glassed_cod))
#    })

    loading_processed_files <- reactive ({
        if (input$ex_btn[1] - btn_counter){
            btn_counter <<- input$ex_btn[1]
            ex <- c("data/abis/eric/3low_freq_fsF.ab1",
                    "data/abis/eric/3low_freq_fsR.ab1")
        }
        calls <- structure("error_reading_Rbin",class = "my_UI_exception")

        if(!is.null(input$goLoadSamples)){
            isolate({
                load_id <- as.numeric(strsplit(input$goLoadSamples, "_")[[1]][2])
            })
        }
        #if(!is.null(input$select_file) || !is.null(ex)) {
        if(!is.null(load_id) || !is.null(ex)) {
            fwd_file      <- g_files[load_id]$FWD_file
            fwd_file_name <- g_files[load_id]$FWD_name
            rev_file      <- g_files[load_id]$REV_file
            rev_file_name <- g_files[load_id]$REV_name
            ref           <- g_files[load_id]$REF
            
            if(!is.null(ex)){
                fwd_file <- fwd_file_name <- ex[1]
                rev_file <- rev_file_name <- ex[2]
                ref = "TP53"
            }
            g_glassed_ref <<- paste("data/refs/",ref,".glassed.intrex.fasta",sep="")
            g_glassed_cod <<- paste("data/refs/",ref,".glassed.codons.rdata",sep="")
            
            if (fwd_file_name == "-"){
                single_rev <- TRUE
                fwd_file <- NULL
            }
            else
                single_rev <- FALSE
            if (rev_file_name == "-")
                rev_file <- NULL
          #  file <- input$select_file$datapath
          #  name <- input$select_file$name
          #  #input$gene_of_interest
          #  full_name <- input$select_file$name
          #  single_rev <- FALSE
#
          #  #getting rid of the date in names
          #  for(i in 1:length(name)){
          #      name[i] <- gsub("_[12][09][0-9][0-9]-[0-1][0-9]-[0123][0-9]_[0-1][0-9]-[0-5][0-9]-[0-5][0-9]","",name[i])
          #      name[i] <- gsub(".abi",".ab1",name[i])
          #  }
          #  #if multiple files uploaded we use the first to
          #  #check if we can distinguish forward and reverse
          #  #otherwise we only take the first file
          #  base <- ""
          #  if(length(name)>=2){
          #      if(gsub("*F.ab1","F",name[1])==gsub("*R.ab1","F",name[2])){
          #          fwd_file <- input$select_file$datapath[1]
          #          fwd_file_name <- full_name[1]
          #          rev_file <- input$select_file$datapath[2]
          #          rev_file_name <- full_name[2]
          #      }else if(gsub("*R.ab1","R",name[1])==gsub("*F.ab1","R",name[2])){
          #          fwd_file <- input$select_file$datapath[2]
          #          fwd_file_name <- full_name[2]
          #          rev_file <- input$select_file$datapath[1]
          #          rev_file_name <- full_name[1]
          #      }else{
          #          fwd_file <- input$select_file$datapath[1]
          #          fwd_file_name <- full_name[1]
          #          rev_file <- NULL
          #          rev_file_name <- "-"
          #      }
          #  }else if(is.null(ex)){ #only one file selected
          #      fwd_file <- input$select_file$datapath
          #      fwd_file_name <- full_name
          #      rev_file <- NULL
          #      rev_file_name <- "-"
          #      ref <- read.fasta(paste0("data/refs/",input$gene_of_interest,".glassed.intrex.fasta"))
          #      re
          #      sm <<- matrix(c(1 ,-1 ,-1 ,-1 ,-1 ,0.5 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,0.1 ,0 ,-1 ,1 ,-1 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.5 ,-1 ,0.1 ,-1 ,0.1 ,0.1 ,0 ,-1 ,-1 ,1 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.5 ,-1 ,0.1 ,0.1 ,-1 ,0.1 ,0 ,-1 ,-1 ,-1 ,1 ,0.5 ,-1 ,-1 ,0.5 ,-1 ,0.5 ,0.1 ,0.1 ,0.1 ,-1 ,0 ,-1 ,-1 ,0.5 ,0.5 ,0.1 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.5 ,0.5 ,-1 ,-1 ,-1 ,0.1 ,0 ,0 ,0 ,0 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.1 ,0.5 ,-1 ,0.5 ,-1 ,0 ,0 ,0.1 ,-1 ,0 ,0 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-1 ,0.5 ,-1 ,0.5 ,0 ,0 ,-1 ,0.1 ,0 ,0 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.5 ,0.5 ,-1 ,0 ,0 ,0 ,0 ,0.1 ,-1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,0.5 ,-1 ,-1 ,0.5 ,0 ,0 ,0 ,0 ,-1 ,0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,-1 ,0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,-0.1 ,0.1 ,0 ,0.1 ,0 ,0 ,0.1 ,0.1 ,0.1 ,-1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0.1 ,0 ,0 ,0.1 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,-1 ,-0.1 ,0.1 ,0.1 ,-0.1 ,0.1 ,-0.1 ,0 ,0 ,0 ,0.1 ,0.1 ,0 ,0 ,0 ,0 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1 ,0.1),15,15,dimnames = list(c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),c("A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N")))
          #      seq<-DNAString(sangerseqR::read.abif(fwd_file)@data$PBAS.2)
          #      pa <- pairwiseAlignment(pattern = ref, subject = seq,type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
          #      score_fwd <- pa@score
          #      pa <- pairwiseAlignment(pattern = ref, subject = reverseComplement(seq),type = "local",substitutionMatrix = sm,gapOpening = -6, gapExtension = -1)
          #      score_rev <- pa@score
          #      if(score_rev>score_fwd){
          #          single_rev <- TRUE
          #      }
          #      #if(!is.null(name)){
          #      #    base <- sapply(strsplit(basename(name),"\\."),
          #      #                   function(x) paste(x[1:(length(x)-1)], collapse="."))
          #      #}
          #  }
            
            #if(substr(base,nchar(base),nchar(base))=="R"){
            if(single_rev){
                isolate({
                    files_info <<- paste0("fwd (F): ",rev_file_name,"\nrev (R): ",fwd_file_name," \n<em>aligned to: ",ref,"</em>",sep="")
                    #single_rev <- TRUE
                })
            }else if(is.null(ex)){
                isolate({
                    files_info <<- paste0("fwd (F): ",fwd_file_name,"\nrev (R): ",rev_file_name," \n<em>aligned to: ",ref,"</em>",sep="")
                })
            }
            if(!is.null(ex)){
                files_info <<- "frameshift deletion (c.277_278delCT) and heterozygous polymorphism (c.215C>G), detectable with these settings:\nmutation minimum peak % =~ 7; minimum quality =~ 16; 'use detected hetero indels' = checked"
                output$files <- renderPrint({cat("<pre>frameshift deletion (c.277_278delCT) and heterozygous polymorphism (c.215C>G), detectable with these settings:\nmutation minimum peak % =~ 7; minimum quality =~ 16; 'use detected hetero indels' = checked</pre>")})
                base = ""
                ex <- NULL
            }

            withProgress(message = paste('processing...',sep=" "), value = 1, {
                # TODO - make sure shiny only allows ab1 files (???)
                #        - otherwise chceck if file has correct format

                tryCatch(
                    g_abif <- sangerseqR::read.abif(fwd_file)@data,
                    error = function(e){output$files <- renderPrint(paste0("<pre>error while reading forward file, are you loading .abi ? ",e$message,"</pre>" ))})
                if(!is.null(rev_file)) {
                    tryCatch(
                        g_abif_rev <- sangerseqR::read.abif(rev_file)@data,
                        error = function(e){output$files <- renderPrint(paste0("<pre>error while reading reverse file, are you loading .abi ? ",e$message,"</pre>" ))})
                }
                else g_abif_rev <- NULL

                res <- NULL
                called <- NULL
                #res <- get_call_data(g_abif,g_abif_rev,input$rm7qual_thres,input$qual_thres,input$aln_min)
                tryCatch(
                    called <- suppressWarnings(get_call_data(g_abif,g_abif_rev,single_rev,g_glassed_ref)),
                    error = function(e){output$files <- renderPrint(paste0("<pre>error while loading calls from abi file : ",e$message,"</pre>" ))})

                if(!is.null(called)){
                    g_minor_het_insertions  <<- NULL
                    g_stored_het_indels     <<- list()
                    g_indels_present        <<- FALSE
                    g_qual_present          <<- called$qual_present
                    intensified  <-  get_intensities(g_abif,g_abif_rev,calls=called$calls,deletions=called$deletions,norm=FALSE,single_rev)
                    g_intens     <<- intensified$intens
                    g_intens_rev <<- intensified$intens_rev
                    calls        <-  annotate_calls(calls=intensified$calls,intens=intensified$intens,intens_rev=intensified$intens_rev,g_glassed_cod)
                    calls        <-  adjust_ref_mut(calls,g_intens_rev)
                    g_max_y      <<- max(c(max(g_intens[,list(A,C,G,T)]),if(is.null(g_intens_rev)) 0 else max(g_intens_rev[,list(A,C,G,T)])))
                    #intrex contains intesities coordinates of start and end of introns/exons with the sequence id (position in sequence coordinates)
                        intrexdat            <- list()
                        intrexdat$intrex     <- list()
                        intrexdat$intrex     <- setnames(calls[!is.na(exon_intron),list(max(id)-min(id)+1,min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","length","trace_peak","end"))
                        intrexdat$intrex     <- setnames(merge(intrexdat$intrex,calls[,list(id,trace_peak)],by="trace_peak"),"trace_peak","start")
                        intrexdat$max_x      <- max(c(nrow(g_intens),nrow(g_intens_rev))) # these numbers should be the same
                        intrexdat$new_sample <- TRUE
                    g_intrexdat       <<- splice_variants(intrexdat)
                    calls             <-  data.table(calls,key="id")
                    g_noisy_neighbors <<- get_noisy_neighbors(calls)
                    if(!called$qual_present){
                        g_files <- paste0(g_files,HTML("\n<strong style=\"color: red;\">no Phred qualities!</strong>"))
                    }
                    files_info <- paste0("<pre>",files_info,"</pre>")
                    output$files      <-  renderPrint({cat(files_info)})
                    g_new_sample      <<- TRUE
                } else return(structure("error_reading_Rbin",class = "my_UI_exception"))
            })
        }
        g_stored_het_indels <<- list()
        g_calls <<- NULL
        return(calls)
    })

    varcall <- reactive({
        if(class(loading_processed_files())[1] != "my_UI_exception") {
            update_chosen_variants()
            goReset_handler()

            calls<-loading_processed_files()
            if(is.null(g_calls)){
                g_calls <<- calls
            }
            #remove added minor het ins
            if(exists("g_minor_het_insertions") && !is.null(g_minor_het_insertions$added)){
                for(i in nrow(g_minor_het_insertions))
                {
                    added <- strsplit(g_minor_het_insertions[i]$added[[1]], split = " ")
                    g_calls <<- g_calls[!(id %in% added[[1]])]
                    ret <- remove_intensities(added,g_calls,g_intens,g_intens_rev,g_intrexdat,g_minor_het_insertions)
                    g_calls      <<- ret$calls
                    g_intens     <<- ret$intens
                    g_intens_rev <<- ret$intens_rev
                    g_intrexdat  <<- ret$intrexdat
                }
                g_minor_het_insertions[,added:=NULL]
                g_minor_het_insertions[,ins_added := NULL]
            }

            g_calls <<- call_variants(g_calls,input$qual_thres_to_call,input$mut_min,input$s2n_min,g_stored_het_indels)
            setkey(g_calls,id)

            report                  <- report_hetero_indels(g_calls)
            g_indels_present       <<- report$indels_present
            g_minor_het_insertions <<- report$minor_het_insertions
            g_hetero_indel_aln     <<- report$hetero_indel_aln
            g_hetero_ins_tab       <<- report$hetero_ins_tab
            g_hetero_del_tab       <<- report$hetero_del_tab
            g_hetero_indel_pid     <<- report$hetero_indel_pid
            g_hetero_indel_report  <<- report$hetero_indel_report
            
            if(input$incorporate_checkbox & g_indels_present){ 
                g_calls <<- incorporate_hetero_indels_func(g_calls,g_hetero_del_tab,g_hetero_ins_tab,g_minor_het_insertions)
            }
            setkey(g_calls,id)

            if(exists("g_minor_het_insertions") && !is.null(g_minor_het_insertions$added)){
                for(i in 1:nrow(g_minor_het_insertions)){
                    ret <- add_intensities(strsplit(g_minor_het_insertions[i]$added[[1]],split= " ")[[1]],g_calls,g_intens,g_intens_rev,g_intrexdat)
                }
                g_calls                          <<- ret$calls
                g_minor_het_insertions$ins_added <<- ret$ins_added
                g_intens                         <<- ret$intens
                g_intens_rev                     <<- ret$intens_rev
                g_intrexdat                      <<- ret$intrexdat  
            }

            g_expected_het_indel <<- get_expected_het_indels(g_calls)
            g_calls   <<- retranslate(g_calls)
            g_choices <<- get_choices(g_calls)
            return(TRUE)

        } else return(FALSE)
    })

    #
    # Render functions reacting to varcall
    #
    output$plot <- renderChromatography({
        if(varcall()) {
            g_intrexdat$max_y <- (g_max_y*100)/input$max_y_p
            ret<-chromatography(g_intens,g_intens_rev,g_intrexdat,g_calls,g_choices,g_new_sample,g_noisy_neighbors,input$show_calls_checkbox,g_qual_present)
            g_new_sample <<- FALSE
            return(ret)
        }
    })

    output$hetero_indel_pid <- renderPrint({
        if(varcall() ) {
            if(!is.null(g_expected_het_indel) && g_expected_het_indel[[1]] * 100 >= 1) het_indel_info <- paste0("expected starting at ~",g_expected_het_indel[[1]] * 100,"% ->\n")
            else                                                                       het_indel_info <- paste0(" ... none expected\n")
            cat(het_indel_info)
            cat(g_hetero_indel_report)
        }
    })

    output$infobox <- renderPrint({
        input$qual_thres_to_call
        #if(loading_processed_files() != "not") {
            if(input$choose_call_pos != "") {
                tryCatch({
                    if(!is.null(g_intens_rev)) {
                        cat(g_calls[id == input$choose_call_pos,paste0("pos ",id,"    ref ",reference,"   user ",user_sample,"    max.peak% ",round(sample_peak_pct,1),"\n ",exon_intron,"  @genomic ",gen_coord,"  @coding ",coding_seq,"  @codon ",codon,"\n\nF  mut ",mut_peak_base_fwd,"  \tQ ",quality_fwd,"  \tpeak% ",round(mut_peak_pct_fwd,digits=1),"  \tS/N ",round(mut_s2n_abs_fwd,digits=1),"\nR  mut ",mut_peak_base_rev,"  \tQ ",quality_rev,"  \tpeak% ",round(mut_peak_pct_rev,digits=1),"  \tS/N ",round(mut_s2n_abs_rev,digits=1),sep="")])
                    } else {
                        cat(g_calls[id == input$choose_call_pos,paste0("pos ",id,"    ref ",reference,"   user ",user_sample,"    max.peak% ",round(sample_peak_pct,1),"\n ",exon_intron,"  @genomic ",gen_coord,"  @coding ",coding_seq,"  @codon ",codon,"\n\n   mut ",mut_peak_base_fwd,"  \tQ ",quality,    "  \tpeak% ",round(mut_peak_pct_fwd,digits=1),"  \tS/N ",round(mut_s2n_abs_fwd,digits=1),sep="")])
                    }
                }, error = function(er){
                    if(grepl("NAs introduced",er)) cat("type an integer number")
                })
                session$sendCustomMessage(type = 'input_change',message = paste0(g_calls[id==input$choose_call_pos]$trace_peak))
            } else cat("")
        #} else cat("load .abi/.ab1 file")  #!todo write this message somehwere else
    })

    #
    # Change single / Change button
    #
    update_chosen_variants <- reactive({
        input$change_btn
        isolate({
            if(!is.null(g_calls[["id"]])){
                if (input$change_user_sample != "") g_calls[id==as.numeric(input$choose_call_pos)]$user_sample <<- input$change_user_sample
                if (input$change_user_mut    != "") g_calls[id==as.numeric(input$choose_call_pos)]$user_mut    <<- input$change_user_mut
                # g_calls[id==as.numeric(input$choose_call_pos)]$user_variant <<- input$change_variant
                g_calls[id==as.numeric(input$choose_call_pos)]$set_by_user <<- TRUE
            }
        })
        return(T)
    })

    #
    # DataTable stuff
    #
#    add_checkboxes <- function(){
#        checkboxes <- paste0('<input type="checkbox" name="row', g_view$id, '" value="', g_view$id, '"',"")
#        for(i in 1:nrow(g_view)) {
#            if(g_view[i]$id %in% g_selected) checkboxes[i] <- paste0(checkboxes[i]," checked>","")
#            else checkboxes[i] <- paste0(checkboxes[i],">","")
#        }
#        return(checkboxes)
#    }

    output$chosen_variants_table <- DT::renderDataTable({
        if(varcall())
        if(!is.null(g_choices) && nrow(g_choices) > 0){
            input$change_btn
            input$reset_btn
            if(varcall() & !is.null(g_choices)) {
                g_view<<-get_view(g_calls,g_choices)
                add_goto_buttons     <- shinyInput(actionButton, g_view$id, 'button_', label = "goto",   onclick = 'Shiny.onInputChange(\"goGoto\",  this.id)' )
                add_reset_buttons    <- shinyInput(actionButton, g_view$id, 'button_', label = "remove", onclick = 'Shiny.onInputChange(\"goReset\",  this.id)' )
                add_lock_buttons     <- shinyInput(actionButton, g_view$id, 'button_', label = NULL,   onclick = 'Shiny.onInputChange(\"goLock\",  this.id)',ico = unlist(lapply(g_view$set_by_user, function(x){if(isTRUE(x)){"lock"}else{ "unlock"}})) )
                cbind(Goto=add_goto_buttons, Reset=add_reset_buttons, Lock=add_lock_buttons, g_view[,list("call position"=id,"genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein,"pri peak %"=sample_peak_pct,"sec peak %"=mut_peak_pct)])
                
            }
        }
    }
#    , options = list(dom = "t",orderClasses=c(-1,-2,-3,-4), paging=F, columnDefs=list(list(targets=c("_all"), searchable=F),list(targets=c(0,1,2,3), orderable=F, title="")))
    , escape=FALSE
    , style= 'bootstrap'
    , options=list("paging"=FALSE,"searching"=FALSE,"ordering"=FALSE,"autoWidth"=FALSE)
    )

#    variant_select <- observe({
#        if(varcall()) {
#            g_selected <<- str_trim(input$rows)
#            g_selected_goto_index <<- 0
#        }
#    })

    #
    # data table buttons
    #
    goGoto_handler <- observe({
        if(is.null(input$goGoto)) return()
        goto_id <- as.numeric(strsplit(input$goGoto, "_")[[1]][2])
        session$sendCustomMessage(type = 'goto',message = paste0(g_calls[id==goto_id]$trace_peak))
        updateTextInput(session,"choose_call_pos",value=paste0(goto_id))
    })

    goReset_handler <- reactive({
        #if(varcall()) {
            if(!is.null(input$goReset)){
                isolate({
                    reset_id <- as.numeric(strsplit(input$goReset, "_")[[1]][2])
                    if(!is.na(g_view[id==reset_id]$ids)){
                        ids <- strsplit(g_view[id==reset_id]$ids," ")[[1]]
                        for( cid in ids){
                            updateTextInput(session,"choose_call_pos",value=paste0(cid))
                            g_calls[id==cid]$user_sample <<- g_calls[id==cid]$reference
                            g_calls[id==cid]$user_mut    <<- g_calls[id==cid]$reference
                            g_calls[id==cid]$set_by_user <<- TRUE
                        }
                    }else{
                        updateTextInput(session,"choose_call_pos",value=paste0(reset_id))
                        #reset button changed to remove variant
                        g_calls[id==reset_id]$user_sample <<- g_calls[id==reset_id]$reference
                        g_calls[id==reset_id]$user_mut    <<- g_calls[id==reset_id]$reference
                        g_calls[id==reset_id]$set_by_user <<- TRUE
                    }
                })
            }
            return(T)
        #}
    })

    goLock_handler <- observe({
        if(is.null(input$goLock)) return()
        isolate({
            lock_id <- as.numeric(strsplit(input$goLock, "_")[[1]][2])
            coding  <- g_view[id==lock_id]$coding
            updateTextInput(session,"choose_call_pos",value=paste0(lock_id))
            if(length(grep("ins|del|dup",coding)) > 0){
                g_stored_het_indels[[coding]] <<- lock_id
            } else {
                g_calls[id==lock_id]$set_by_user <<- TRUE
            }
            #output$goLock <- renderUI({actionButton(input$goLock, icon = icon("lock"))})

        })
    })

    #
    # Input binding from JS
    #
    goClick_handler <- observe({
        if(is.null(input$pos_click)) return()
        isolate({
            updateTextInput(session,"choose_call_pos",value=paste0(input$pos_click$id))
        })
    })

    #
    # Send message to JS
    #
    s2n_slider_handler <- observe({
            session$sendCustomMessage(type = 's2n_min',message = paste0(input$s2n_min))
    })


    join_traces <- observe({
        loading_processed_files()
        if(is.null(g_intens_rev))   session$sendCustomMessage(type = "join",message = "TRUE")
        else                        session$sendCustomMessage(type = "join",message = paste0(input$join_traces_checkbox))
    })

#     show_calls <- observe({
#         if(varcall()){
#             session$sendCustomMessage(type = "show",message = paste0(input$show_calls_checkbox))
#         }
#     })

    set_opacity <- observe({
        if(varcall()){
            opac_fwd <- 1 + (input$opacity/100)
            opac_rev <- 1 - (input$opacity/100)
            if(opac_fwd>1)opac_fwd <- 1
            if(opac_rev>1)opac_rev <- 1
            session$sendCustomMessage(type = "opac_f",message = paste0(opac_fwd))
            session$sendCustomMessage(type = "opac_r",message = paste0(opac_rev))
        }
    })


    #
    # EXPORT
    #
    output$export_btn <- downloadHandler(
        filename = function() {
            paste('data-', Sys.Date(), '.xlsx', sep='')
        },
        content = function(con) {
            out<-data.table("genomic coordinate"=character(),"coding variant"=character(),"protein variant"=character())
            g_selected <- g_view$id[as.numeric(input$chosen_variants_table_rows_selected)]
            for(i in 1:nrow(g_view)) {
                if(g_view[i]$id %in% g_selected) out<-rbind(out,g_view[i,list("genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein)])
            }
            if(length(g_selected)==0)
                write.xlsx(g_view[,list("genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein)], con)
            else write.xlsx(out,con)
            
        }
    )

    #
    # Other tabs
    #

#    output$aln <- renderPrint({
#      if(varcall() ) {
#          cat("(P)rimary vs (S)econdary (or consensus of fwd+rev secondaries)\n\nidentified insertions:\n")
#          if(is.na(g_hetero_ins_tab[1])) cat("no insertions\n")
#          else print(g_hetero_ins_tab)
#          cat("\nidentified deletions:\n")
#          if(is.na(g_hetero_del_tab[1])) cat("no deletions\n")
#          else print(g_hetero_del_tab)
#          cat("\n")
#          writePairwiseAlignments(g_hetero_indel_aln, block.width = 150)
#      }
#    })
#
#    output$het_histogram <- renderPlot({
#        if(varcall() ) {
#            if(!is.null(g_expected_het_indel)){
#                plot(g_expected_het_indel$hist)
#                if(g_expected_het_indel$min != 0) abline(v = g_expected_het_indel$min,col = "red")
#                #if(g_expected_het_indel$max != 0) abline(v = g_expected_het_indel$max,col = "orange")
#            }
#        }
#    },height = 600)
#
#     output$call_table <- shiny::renderDataTable({
#         if(varcall() & !is.null(g_calls)) { g_calls }
#     }
#     ,options = list(paging=F
#                     ,columnDefs=list(list(searchable=F, orderable=F, title=""))
#                     ))
#
#     output$intens_table <- shiny::renderDataTable({
#         if(varcall() & !is.null(g_intens)) { g_intens }
#     }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))

#     output$intens_table_rev <- shiny::renderDataTable({
#         if(varcall() & !is.null(g_intens_rev)) { g_intens_rev }
#     }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    #used for conditional display
    output$reverse <- reactive({
        loading_processed_files()
        return(!is.null(g_intens_rev))
    })

    output$indels_present <- reactive({
        varcall()
        return(g_indels_present)
    })
    #output$show_sample_brows <- reactive({
    #    input$mng_samples_btn
    #    return(FALSE)
    #})
    #not shure why but I need this here to make it work (regirsters the output variable?)
    outputOptions(output, 'reverse', suspendWhenHidden=FALSE)
    outputOptions(output, 'indels_present', suspendWhenHidden=FALSE)
    #outputOptions(output, 'show_sample_brows',suspendWhenHidden=FALSE)
})
