# library(shiny)
library(sangerseqR)
source("procAbi.R")
source("helpers.R")

g_calls                 <<- NULL             #annotated basecall data
makeReactiveBinding("g_calls")
g_intens                <<- NULL             #intensities file
g_intens_rev            <<- NULL             #optional reverse intensities file
g_helperdat             <<- NULL             #helper data used in graphs
g_choices               <<- NULL
g_selected              <<- NULL
g_selected_goto_index   <<- 0
g_max_y                 <<- NULL
g_varcall               <<- FALSE
g_hetero_calls          <<- 0
g_hetero_indel_pid      <<- 0
g_hetero_ins_tab        <<- NULL
g_hetero_del_tab        <<- NULL

shinyServer(function(input,output,session) {

    get_file <- reactive({
       # if (!is.null(input$select_file)) return(input$select_file$datapath)
        if (!is.null(input$select_file)) return(input$select_file)
        else return(NULL)
    })

    loading_processed_files <- reactive ({

        if(!is.null(get_file())) {
            ret <- "not"
            g_varcall <<- FALSE
            file <- get_file()$datapath
            name <- get_file()$name

            #if multiple files uploaded we use the first to
            #check if we can distinguish forward and reverse
            #otherwise we only take the first file
            if(length(name)>=2){
                if(gsub("F.*","F",name[1])==gsub("R.*","F",name[2])){
                    fwd_file <- get_file()$datapath[1]
                    fwd_file_name <- name[1]
                    rev_file <- get_file()$datapath[2]
                    rev_file_name <- name[2]
                }else if(gsub("R.*","R",name[1])==gsub("F.*","R",name[2])){
                    fwd_file <- get_file()$datapath[2]
                    fwd_file_name <- name[2]
                    rev_file <- get_file()$datapath[1]
                    rev_file_name <- name[1]
                }else{
                    fwd_file <- get_file()$datapath[1]
                    fwd_file_name <- name[1]
                    rev_file <- NULL
                    rev_file_name <- "-"
                }
            }else{
                fwd_file <- get_file()$datapath
                fwd_file_name <- name
                rev_file <- NULL
                rev_file_name <- "-"
            }
            g_files <<- paste0("fwd: ",fwd_file_name,"\nrev: ",rev_file_name,sep="")

            withProgress(message = paste('processing...',sep=" "), value = 1, {
                # TODO - make sure shiny only allows ab1 files (???)
                #        - otherwise chceck if file has correct format

                tryCatch(
                    g_abif <- sangerseqR::read.abif(fwd_file)@data,
                    error = function(e){output$infobox <- renderPrint(paste0("error while reading forward file, are you loading .abi ? ",e$message ))})
                if(!is.null(rev_file)) {
                    tryCatch(
                        g_abif_rev <- sangerseqR::read.abif(rev_file)@data,
                        error = function(e){output$infobox <- renderPrint(paste0("error while reading reverse file, are you loading .abi ? ",e$message ))})
                }
                else g_abif_rev <- NULL

                res <- NULL
                called <- NULL
                #res <- get_call_data(g_abif,g_abif_rev,input$rm7qual_thres,input$qual_thres,input$aln_min)
                tryCatch(
                    called <- suppressWarnings(get_call_data(g_abif,g_abif_rev)),
                    error = function(e){output$infobox <- renderPrint(paste0("error while loading calls from abi file : ",e$message ))})

                if(!is.null(called)){
                    intensified         <-  get_intensities(g_abif,g_abif_rev,calls=called$calls,deletions=called$deletions,norm=FALSE)
                    calls               <-  annotate_calls(calls=intensified$calls,intens=intensified$intens,intens_rev=intensified$intens_rev)
                    g_intens            <<- intensified$intens
                    g_intens_rev        <<- intensified$intens_rev
                    g_max_y             <<- max(c(max(g_intens[,list(A,C,G,T)]),if(is.null(g_intens_rev)) 0 else max(g_intens_rev[,list(A,C,G,T)])))
                    #helper_intrex contains intesities coordinates of start and end of introns/exons with the sequence id (position in sequence coordinates)
                        helperdat               <- list()
                        helperdat$helper_intrex <- list()
                        helperdat$helper_intrex <- setnames(calls[!is.na(exon_intron),list(min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","trace_peak","end"))
                        helperdat$helper_intrex <- setnames(merge(helperdat$helper_intrex,calls[,list(id,trace_peak)],by="trace_peak"),"trace_peak","start")
                        helperdat$max_x         <- max(c(nrow(g_intens),nrow(g_intens_rev))) #although these numbers should be the same
                        helperdat$new_sample    <- TRUE
                    g_helperdat         <<- helperdat
                    g_calls             <<- data.table(calls,key="id")
                    ret<-"loaded"
                    g_new_sample <<- TRUE
                }
            })
            return(ret)
        }
        return("not")
    })
    output$files <- renderPrint({
        if(loading_processed_files() != "not") {
            cat(g_files)
        }
    })

    varcall <- reactive({
        if(loading_processed_files() != "not"){
            g_calls <<- call_variants(g_calls,input$qual_thres_to_call,input$mut_min,input$s2n_min)

            rep <- report_hetero_indels(g_calls)
            g_hetero_indel_aln <<- rep[[1]]
            g_hetero_indel_pid <<- rep[[2]]
            g_hetero_ins_tab   <<- rep[[3]]
            g_hetero_del_tab   <<- rep[[4]]

            if(input$incorporate_checkbox) g_calls <<- incorporate_hetero_indels_func(g_calls)

            g_calls   <<- retranslate(g_calls)
            g_choices <<- get_choices(g_calls)

            g_varcall <<- TRUE
        }
        return(g_varcall)
    })

    output$plot <- renderChromatography({
        if((loading_processed_files() != "not") & varcall() ) {
            g_helperdat$max_y <- (g_max_y*100)/input$max_y_p
            ret<-chromatography(g_intens,g_intens_rev,g_helperdat,g_calls,g_choices,g_new_sample)
            g_new_sample <<- FALSE
            return(ret)
        }
    })

    output$aln <- renderPrint({
        if(loading_processed_files() != "not" & varcall() ) {
            cat("primary vs secondary (or consensus of fwd+rev secondaries)\n\nidentified insertions:\n")
            if(is.na(g_hetero_ins_tab[1])) cat("no insertions\n")
            else print(g_hetero_ins_tab)
            cat("\nidentified deletions:\n")
            if(is.na(g_hetero_del_tab[1])) cat("no deletions\n")
            else print(g_hetero_del_tab)
            cat("\n")
            writePairwiseAlignments(g_hetero_indel_aln, block.width = 150)
        }
    })
    output$hetero_indel_pid <- renderPrint({
        if(loading_processed_files() != "not" & varcall() ) {
            cat(g_hetero_indel_pid)
        }
    })
    output$hetero_indel_tab <- renderPrint({
        if(loading_processed_files() != "not" & varcall() ) {
            cat((g_hetero_ins_tab[,2]-g_hetero_ins_tab[,1]+1),"/",(g_hetero_del_tab[,2]-g_hetero_del_tab[,1]+1))
        }
    })

#     incorporate_hetero_indels <- observe({
#         input$incorporate_btn
#             if(isolate(loading_processed_files()) != "not") {
#                 g_calls <<- incorporate_hetero_indels_func(g_calls)
#             }
#     })

    output$infobox <- renderPrint({
        if(loading_processed_files() != "not") {
            if(input$choose_call_pos != "") {
                tryCatch({
                    if(!is.null(g_intens_rev)) {
                        cat(g_calls[id == input$choose_call_pos,paste0("pos ",id,"   ref ",reference,"   call ",cons,"   user ",user_sample,"  max.peak% ",round(sample_peak_pct,1),"\n",exon_intron,"  @genomic ",gen_coord,"  @coding ",coding_seq,"  @codon ",codon,"\nfwd mut ",mut_peak_base_fwd,"  \tpeak% ",round(mut_peak_pct_fwd,digits=1),"  \tS/N ",round(mut_s2n_abs_fwd,digits=1),"\nrev mut ",mut_peak_base_rev,"  \tpeak% ",round(mut_peak_pct_rev,digits=1),"  \tS/N ",round(mut_s2n_abs_rev,digits=1),sep="")])
                    } else {
                        cat(g_calls[id == input$choose_call_pos,paste0("pos ",id,"   ref ",reference,"   call ",cons,"   user ",user_sample,"  max.peak% ",round(sample_peak_pct,1),"\n",exon_intron,"  @genomic ",gen_coord,"  @coding ",coding_seq,"  @codon ",codon,"\nfwd mut ",mut_peak_base_fwd,"  \tpeak% ",round(mut_peak_pct_fwd,digits=1),"  \tS/N ",round(mut_s2n_abs_fwd,digits=1),sep="")])
                    }
                }, error = function(er){
                    if(grepl("NAs introduced",er)) cat("type an integer number")
                })
            } else cat("")
        } else cat("load .abi/.ab1 file")
    })

    update_chosen_variants <- observe({
        input$change_btn
        isolate({
            if(loading_processed_files() != "not") {
                if (input$change_user_sample != "") g_calls[id==as.numeric(input$choose_call_pos)]$user_sample <<- input$change_user_sample
                if (input$change_user_mut    != "") g_calls[id==as.numeric(input$choose_call_pos)]$user_mut    <<- input$change_user_mut
                # g_calls[id==as.numeric(input$choose_call_pos)]$user_variant <<- input$change_variant
                g_calls[id==as.numeric(input$choose_call_pos)]$set_by_user <<- TRUE
            }
        })
    })

    add_checkboxes <- function(){
        checkboxes <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '"',"")
        for(i in 1:nrow(g_choices)) {
            if(g_choices[i]$id %in% g_selected) checkboxes[i] <- paste0(checkboxes[i]," checked>","")
            else checkboxes[i] <- paste0(checkboxes[i],">","")
        }
        return(checkboxes)
    }

    output$chosen_variants_table <- shiny::renderDataTable({
        if(varcall())
        if(!is.null(g_choices) && nrow(g_choices) > 0){
            input$change_btn
            input$reset_btn
            if(loading_processed_files() != "not" & !is.null(g_choices)) {

                #add_checkbox_buttons <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '">',"")
                #add_edit_buttons <- paste0('<a class="go-edit" href="" data-id="', g_choices$id, '"><i class="fa fa-crosshairs"></i></a>')
                add_checkbox_buttons <- add_checkboxes()
                add_goto_buttons     <- paste0('<input type="button" class="go-goto"  value="goto"  name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
                add_reset_buttons    <- paste0('<input type="button" class="go-reset" value="reset" name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
                add_lock_buttons     <- paste0('<input type="button" class="go-lock"  value="lock"  name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")

                # cbind(Pick=add_checkbox_buttons, Edit=add_edit_buttons, Zoom=add_zoom_buttons, g_choices[,list("call position"=id,"coding variant"=coding,"protein variant"=protein,reference,"sample variant"=user_sample,"%"=sample_peak_pct,"mutant variant"=user_mut,"%"=mut_peak_pct)])
                cbind(Pick=add_checkbox_buttons, Goto=add_goto_buttons, Reset=add_reset_buttons, Lock=add_lock_buttons, g_choices[,list("call position"=id,"genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein)])

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

    output$call_table <- shiny::renderDataTable({
    	if(loading_processed_files() != "not" & !is.null(g_calls)) { g_calls }
    }
    ,options = list(paging=F, columnDefs=list(list(searchable=F, orderable=F, title=""))))

#     output$intens_table <- shiny::renderDataTable({
#         if(loading_processed_files() != "not" & !is.null(g_intens)) { g_intens }
#     }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))
#
#     output$intens_table_rev <- shiny::renderDataTable({
#         if(loading_processed_files() != "not" & !is.null(g_intens_rev)) { g_intens_rev }
#     }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    variant_select <- observe({
        if(loading_processed_files() != "not") {
            g_selected <<- str_trim(input$rows)
            g_selected_goto_index <<- 0
        }
    })

    goGoto_handler <- observe({
        if(loading_processed_files() != "not") {
            if(is.null(input$goGoto)) return()
            session$sendCustomMessage(type = 'goto',message = paste0(g_calls[id==input$goGoto$id]$trace_peak))
            updateTextInput(session,"choose_call_pos",value=paste0(input$goGoto$id))
        }
    })
    goReset_handler <- observe({
        if(loading_processed_files() != "not") {
            if(is.null(input$goReset)) return()
            isolate({
                updateTextInput(session,"choose_call_pos",value=paste0(input$goReset$id))
                g_calls[id==input$goReset$id]$user_sample <<- g_calls[id==input$goReset$id]$cons
                g_calls[id==input$goReset$id]$user_mut    <<- g_calls[id==input$goReset$id]$cons
                g_calls[id==input$goReset$id]$set_by_user <<- TRUE
            })
        }
    })
    goLock_handler <- observe({
        if(loading_processed_files() != "not") {
            if(is.null(input$goLock)) return()
            isolate({
                updateTextInput(session,"choose_call_pos",value=paste0(input$goLock$id))
                g_calls[id==input$goLock$id]$set_by_user <<- TRUE
            })
        }
    })
    goClick_handler <- observe({
        if(loading_processed_files() != "not") {
            if(is.null(input$posClick)) return()
            updateTextInput(session,"choose_call_pos",value=paste0(input$posClick$id))
        }
    })

    s2n_slider_handler <- observe({
        if(loading_processed_files()!= "not"){
            session$sendCustomMessage(type = 'mut_min',message = paste0(input$s2n_min))
        }
    })

#     reset_handler <- observe({
#         input$reset_btn
#         isolate({
#             if(length(g_selected) != 0) {
#                 #g_choices <<- g_choices[-match(as.numeric(g_selected),gid)]
#                 for(i in as.numeric(g_selected)) {
#                     g_calls[id==i,]$user_sample <<- g_calls[id==i,]$reference
#                     g_calls[id==i,]$set_by_user <<- TRUE
#                 }
#                 #g_selected <<-  NULL
#             }
#         })
#     })

    split_traces <- observe({
        if(loading_processed_files() != "not"){
            session$sendCustomMessage(type = "split",message = paste0(input$offset_traces_checkbox))
        }
    })
    set_opacity <- observe({
        if(loading_processed_files() != "not"){
            opac_fwd <- 1 + (input$opacity/100)
            opac_rev <- 1 - (input$opacity/100)
            #if(opac_fwd>1)opac_fwd <- 1
            #if(opac_rev>1)opac_rev <- 1
            session$sendCustomMessage(type = "opac_f",message = paste0(opac_fwd))
            session$sendCustomMessage(type = "opac_r",message = paste0(opac_rev))
        }
    })
})
