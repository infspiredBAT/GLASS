# library(shiny)
library(sangerseqR)
library(xlsx)
source("procAbi.R")
source("helpers.R")

#options(shiny.reactlog=TRUE)
g_calls                 <<- NULL             #annotated basecall data
#makeReactiveBinding("g_calls")
g_intens                <<- NULL             #intensities file
g_intens_rev            <<- NULL             #optional reverse intensities file
g_intrexdat             <<- NULL             #intrex data used in graphs
g_choices               <<- NULL
g_noisy_neighbors       <<- NULL
g_view                  <<- NULL
g_selected              <<- NULL
g_selected_goto_index   <<- 0
g_max_y                 <<- NULL
g_hetero_calls          <<- 0
g_hetero_indel_pid      <<- 0
g_hetero_ins_tab        <<- NULL
g_hetero_del_tab        <<- NULL

shinyServer(function(input,output,session) {

#     get_file <- reactive({
#        # if (!is.null(input$select_file)) return(input$select_file$datapath)
#         if (!is.null(input$select_file)) return(input$select_file)
#         else return(NULL)
#     })

    loading_processed_files <- reactive ({

        calls <- structure("error_reading_Rbin",class = "my_UI_exception")
        if(!is.null(input$select_file)) {
            file <- input$select_file$datapath
            name <- input$select_file$name
            single_rev <- FALSE

            #if multiple files uploaded we use the first to
            #check if we can distinguish forward and reverse
            #otherwise we only take the first file
            if(length(name)>=2){
                if(gsub("F.*","F",name[1])==gsub("R.*","F",name[2])){
                    fwd_file <- input$select_file$datapath[1]
                    fwd_file_name <- name[1]
                    rev_file <- input$select_file$datapath[2]
                    rev_file_name <- name[2]
                }else if(gsub("R.*","R",name[1])==gsub("F.*","R",name[2])){
                    fwd_file <- input$select_file$datapath[2]
                    fwd_file_name <- name[2]
                    rev_file <- input$select_file$datapath[1]
                    rev_file_name <- name[1]
                }else{
                    fwd_file <- input$select_file$datapath[1]
                    fwd_file_name <- name[1]
                    rev_file <- NULL
                    rev_file_name <- "-"
                }
            }else{ #only one file selected
                fwd_file <- input$select_file$datapath
                fwd_file_name <- name
                rev_file <- NULL
                rev_file_name <- "-"
            }
            base <- sapply(strsplit(basename(fwd_file_name),"\\."),
                           function(x) paste(x[1:(length(x)-1)], collapse="."))
            if(substr(base,nchar(base),nchar(base))=="R"){
                g_files <<- paste0("fwd (F): ",rev_file_name,"\nrev (R): ",fwd_file_name,sep="")
                single_rev <- TRUE
            }else{
                g_files <<- paste0("fwd (F): ",fwd_file_name,"\nrev (R): ",rev_file_name,sep="")
            }

            withProgress(message = paste('processing...',sep=" "), value = 1, {
                # TODO - make sure shiny only allows ab1 files (???)
                #        - otherwise chceck if file has correct format

                tryCatch(
                    g_abif <- sangerseqR::read.abif(fwd_file)@data,
                    error = function(e){output$files <- renderPrint(paste0("error while reading forward file, are you loading .abi ? ",e$message ))})
                if(!is.null(rev_file)) {
                    tryCatch(
                        g_abif_rev <- sangerseqR::read.abif(rev_file)@data,
                        error = function(e){output$files <- renderPrint(paste0("error while reading reverse file, are you loading .abi ? ",e$message ))})
                }
                else g_abif_rev <- NULL

                res <- NULL
                called <- NULL
                #res <- get_call_data(g_abif,g_abif_rev,input$rm7qual_thres,input$qual_thres,input$aln_min)
                glassed_ref <- paste("data/refs/",input$gene_of_interest,".glassed.intrex.fasta",sep="")
                glassed_cod <- paste("data/refs/",input$gene_of_interest,".glassed.codons.rdata",sep="")
                tryCatch(
                    called <- suppressWarnings(get_call_data(g_abif,g_abif_rev,single_rev,glassed_ref)),
                    error = function(e){output$files <- renderPrint(paste0("error while loading calls from abi file : ",e$message ))})

                if(!is.null(called)){
                    intensified  <-  get_intensities(g_abif,g_abif_rev,calls=called$calls,deletions=called$deletions,norm=FALSE,single_rev)
                    calls        <-  annotate_calls(calls=intensified$calls,intens=intensified$intens,intens_rev=intensified$intens_rev,glassed_cod)
                    g_intens     <<- intensified$intens
                    g_intens_rev <<- intensified$intens_rev
                    g_max_y      <<- max(c(max(g_intens[,list(A,C,G,T)]),if(is.null(g_intens_rev)) 0 else max(g_intens_rev[,list(A,C,G,T)])))
                    #intrex contains intesities coordinates of start and end of introns/exons with the sequence id (position in sequence coordinates)
                        intrexdat            <- list()
                        intrexdat$intrex     <- list()
                        intrexdat$intrex     <- setnames(calls[!is.na(exon_intron),list(max(id)-min(id)+1,min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","length","trace_peak","end"))
                        intrexdat$intrex     <- setnames(merge(intrexdat$intrex,calls[,list(id,trace_peak)],by="trace_peak"),"trace_peak","start")
                        intrexdat$max_x      <- max(c(nrow(g_intens),nrow(g_intens_rev))) #although these numbers should be the same
                        intrexdat$new_sample <- TRUE
                    g_intrexdat       <<- splice_variants(intrexdat)
                    calls             <-  data.table(calls,key="id")
                    g_noisy_neighbors <<- get_noisy_neighbors(calls)
                    output$files      <-  renderPrint({cat(g_files)})
                    g_new_sample      <<- TRUE
                } else return(structure("error_reading_Rbin",class = "my_UI_exception"))
            })
        }
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
                setkey(g_calls,id)
            }
            g_calls <<- call_variants(g_calls,input$qual_thres_to_call,input$mut_min,input$s2n_min)

            rep <- report_hetero_indels(g_calls)
            g_hetero_indel_aln <<- rep[[1]]
            g_hetero_indel_pid <<- rep[[2]]
            g_hetero_ins_tab   <<- rep[[3]]
            g_hetero_del_tab   <<- rep[[4]]
            if(input$incorporate_checkbox) g_calls <<- incorporate_hetero_indels_func(g_calls)

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
            ret<-chromatography(g_intens,g_intens_rev,g_intrexdat,g_calls,g_choices,g_new_sample,g_noisy_neighbors,input$show_calls_checkbox)
            g_new_sample <<- FALSE
            return(ret)
        }
    })

    output$hetero_indel_pid <- renderPrint({
        if(varcall() ) {
            cat(g_hetero_indel_pid,"%\n",(g_hetero_ins_tab[,2]-g_hetero_ins_tab[,1]+1),"/",(g_hetero_del_tab[,2]-g_hetero_del_tab[,1]+1))
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
    add_checkboxes <- function(){
        checkboxes <- paste0('<input type="checkbox" name="row', g_view$id, '" value="', g_view$id, '"',"")
        for(i in 1:nrow(g_view)) {
            if(g_view[i]$id %in% g_selected) checkboxes[i] <- paste0(checkboxes[i]," checked>","")
            else checkboxes[i] <- paste0(checkboxes[i],">","")
        }
        return(checkboxes)
    }

    output$chosen_variants_table <- shiny::renderDataTable({
        if(varcall())
        if(!is.null(g_choices) && nrow(g_choices) > 0){
            input$change_btn
            input$reset_btn
            if(varcall() & !is.null(g_choices)) {

                g_view<<-get_view(g_choices)
                #add_checkbox_buttons <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '">',"")
                #add_edit_buttons <- paste0('<a class="go-edit" href="" data-id="', g_choices$id, '"><i class="fa fa-crosshairs"></i></a>')
                add_checkbox_buttons <- add_checkboxes()
                add_goto_buttons     <- paste0('<input type="button" class="go-goto"  value="goto"  name="btn',g_view$id,'" data-id="',g_view$id,'"',">")
                add_reset_buttons    <- paste0('<input type="button" class="go-reset" value="remove" name="btn',g_view$id,'" data-id="',g_view$id,'"',">")
                add_lock_buttons     <- paste0('<input type="button" class="go-lock"  value="lock"  name="btn',g_view$id,'" data-id="',g_view$id,'"',">")

                # cbind(Pick=add_checkbox_buttons, Edit=add_edit_buttons, Zoom=add_zoom_buttons, g_choices[,list("call position"=id,"coding variant"=coding,"protein variant"=protein,reference,"sample variant"=user_sample,"%"=sample_peak_pct,"mutant variant"=user_mut,"%"=mut_peak_pct)])

                cbind(Pick=add_checkbox_buttons, Goto=add_goto_buttons, Reset=add_reset_buttons, Lock=add_lock_buttons, g_view[,list("call position"=id,"genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein,"pri peak %"=sample_peak_pct,"sec peak %"=mut_peak_pct)])

            }
        }
    }
    , options = list(dom = "t",orderClasses=c(-1,-2,-3,-4), paging=F, columnDefs=list(list(targets=c("_all"), searchable=F),list(targets=c(0,1,2,3), orderable=F, title="")))
    , escape=c(-1,-2,-3,-4)
    , callback =
        "function(table) {
            table.on('change.dt', 'tr td input:checkbox', function() {
                setTimeout(function () {
                    Shiny.onInputChange('rows', $(this).add('tr td input:checkbox:checked').parent().siblings(':nth-child(4)').map(function() {
                        return $(this).text();
                    }).get())
                }, 10);
            });
            table.on('click', '.go-goto', function(e) {
                e.preventDefault();
                $el = $(this);
                var id_data = $el.data('id');
                Shiny.onInputChange('goGoto', {
                    id: id_data
                });
            });
            table.on('click', '.go-reset', function(e) {
                e.preventDefault();
                $el = $(this);
                var id_data = $el.data('id');
                Shiny.onInputChange('goReset', {
                    id: id_data
                });
            });
            table.on('click', '.go-lock', function(e) {
                e.preventDefault();
                $el = $(this);
                var id_data = $el.data('id');
                Shiny.onInputChange('goLock', {
                    id: id_data
                });
            });
        }"
    )

    variant_select <- observe({
        if(varcall()) {
            g_selected <<- str_trim(input$rows)
            g_selected_goto_index <<- 0
        }
    })

    #
    # data table buttons
    #
    goGoto_handler <- observe({
        if(is.null(input$goGoto)) return()
        session$sendCustomMessage(type = 'goto',message = paste0(g_calls[id==input$goGoto$id]$trace_peak))
        updateTextInput(session,"choose_call_pos",value=paste0(input$goGoto$id))
    })

    goReset_handler <- reactive({
        #if(varcall()) {
            if(!is.null(input$goReset)){
                isolate({
                    if(!is.na(g_view[id==input$goReset$id]$ids)){
                        ids <- strsplit(g_view[id==input$goReset$id]$ids," ")[[1]]
                        for( cid in ids){
                            updateTextInput(session,"choose_call_pos",value=paste0(cid))
                            g_calls[id==cid]$user_sample <<- g_calls[id==cid]$reference
                            g_calls[id==cid]$user_mut    <<- g_calls[id==cid]$reference
                            g_calls[id==cid]$set_by_user <<- TRUE
                        }
                    }else{
                        updateTextInput(session,"choose_call_pos",value=paste0(input$goReset$id))
                        #reset button changed to remove variant
                        g_calls[id==input$goReset$id]$user_sample <<- g_calls[id==input$goReset$id]$reference
                        g_calls[id==input$goReset$id]$user_mut    <<- g_calls[id==input$goReset$id]$reference
                        g_calls[id==input$goReset$id]$set_by_user <<- TRUE
                    }
                })
            }
            return(T)
        #}
    })

    goLock_handler <- observe({
        if(is.null(input$goLock)) return()
        isolate({
            updateTextInput(session,"choose_call_pos",value=paste0(input$goLock$id))
            g_calls[id==input$goLock$id]$set_by_user <<- TRUE
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
            #if(opac_fwd>1)opac_fwd <- 1
            #if(opac_rev>1)opac_rev <- 1
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
            out<-data.table()
            for(i in 1:nrow(g_view)) {
                if(g_view[i]$id %in% g_selected) out<-rbind(out,g_view[id==i,list("genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein)])
            }
            write.xlsx(g_view[,list("genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein)], con)
        }
    )

    #
    # Other tabs
    #

    output$aln <- renderPrint({
      if(varcall() ) {
          cat("(P)rimary vs (S)econdary (or consensus of fwd+rev secondaries)\n\nidentified insertions:\n")
          if(is.na(g_hetero_ins_tab[1])) cat("no insertions\n")
          else print(g_hetero_ins_tab)
          cat("\nidentified deletions:\n")
          if(is.na(g_hetero_del_tab[1])) cat("no deletions\n")
          else print(g_hetero_del_tab)
          cat("\n")
          writePairwiseAlignments(g_hetero_indel_aln, block.width = 150)
      }
    })

    output$call_table <- shiny::renderDataTable({
        if(varcall() & !is.null(g_calls)) { g_calls }
    }
    ,options = list(paging=F, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    output$intens_table <- shiny::renderDataTable({
        if(varcall() & !is.null(g_intens)) { g_intens }
    }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    output$intens_table_rev <- shiny::renderDataTable({
        if(varcall() & !is.null(g_intens_rev)) { g_intens_rev }
    }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    #used for conditional display of options for double strand data
    output$reverse <- reactive({
        loading_processed_files()
        return(!is.null(g_intens_rev))
    })
    #not shure why but I need this here to make it work (regirsters the output variable?)
    outputOptions(output, 'reverse', suspendWhenHidden=FALSE)
})
