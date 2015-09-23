# library(shiny)
library(sangerseqR)
source("procAbi.R")
source("helpers.R")

g_calls         <- NULL             #annotated basecall data
makeReactiveBinding("g_calls")
g_intens        <- NULL             #intensities file
g_intens_rev    <- NULL             #optional reverse intensities file
g_helperdat     <- NULL             #helper data used in graphs
g_choices       <- NULL
g_selected      <- NULL
g_selected_zoom_index <- 0
g_max_y         <- NULL
g_varcall       <- FALSE

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
                    rev_file <- get_file()$datapath[2]
                }else if(gsub("R.*","R",name[1])==gsub("F.*","R",name[2])){
                    fwd_file <- get_file()$datapath[2]
                    rev_file <- get_file()$datapath[1]
                }else{
                    fwd_file <- get_file()$datapath[1]
                    rev_file <- NULL
                }
            }else{
                fwd_file <- get_file()$datapath
                rev_file <- NULL
            }

            withProgress(message = paste('processing...',sep=" "), value = 1, {
                # TODO - make sure shiny only allows ab1 files (???)
                #        - otherwise chceck if file has correct format
                
                tryCatch(
                    g_abif <- sangerseqR::read.abif(fwd_file)@data,
                    error = function(e){output$infobox <- renderPrint(paste0("An error occured while reading forward input file, make sure you load an .abi file: ",e$message ))})
                if(!is.null(rev_file)) {
                    tryCatch(
                        g_abif_rev <- sangerseqR::read.abif(rev_file)@data,
                        error = function(e){output$infobox <- renderPrint(paste0("An error occured while reading reverse input file, make sure you load an .abi file: ",e$message ))})
                }
                else g_abif_rev <- NULL
                
                res <- NULL
                #res <- get_call_data(g_abif,g_abif_rev,input$rm7qual_thres,input$qual_thres,input$aln_min)
                tryCatch(
                    called <- suppressWarnings(get_call_data(g_abif,g_abif_rev)),
                    error = function(e){output$infobox <- renderPrint(paste0("An error occured while loading calls from abi file with the following error message: ",e$message ))})

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

    varcall <- reactive({
        if(loading_processed_files() != "not"){
            g_calls     <<- call_variants(g_calls,input$qual_thres,input$mut_min,input$s2n_min)
            g_calls     <<- retranslate(g_calls)
            
            g_choices   <<- g_calls[user_sample != reference & trace_peak != "NA" & !is.na(gen_coord)]
            g_varcall   <<- TRUE
        }
        return(g_varcall)
    })

    output$plot <- renderChromatography({
        #lg_calls
        if((loading_processed_files() != "not") & varcall() ) {
            #g_choices <<- call_variants
            g_helperdat$max_y <- (g_max_y*100)/input$max_y_p
            ret<-chromatography(g_intens,g_intens_rev,g_helperdat,g_calls,g_choices,g_new_sample)
            g_new_sample <<- FALSE
            return(ret)
        }
    })

    output$aln <- renderPrint({
        if(loading_processed_files() != "not" & varcall() ) {
            return(report_hetero_indels(g_calls))
        }
    })

#     output$muts <- renderPrint({
#         if(loading_processed_files() != "not" & varcall() ) {
#             cat("soon")
#         }
#     })

    output$infobox <- renderPrint({
        if(loading_processed_files() != "not") {
            if(input$choose_variance != "") {
                tryCatch({
                    if(!is.null(g_intens_rev)) {
                        cat(g_calls[id == input$choose_variance,paste0("ref ",reference,"   user ",user_sample,"   orig ",cons,"\n",exon_intron,"   @ ",gen_coord,"\nfwd mut ",mut_peak_base_fwd,"  \tpeak% ",round(mut_peak_pct_fwd,digits=1),"  \tS/N ",round(mut_s2n_abs_fwd,digits=1),"\nrev mut ",mut_peak_base_rev,"  \tpeak% ",round(mut_peak_pct_rev,digits=1),"  \tS/N ",round(mut_s2n_abs_rev,digits=1),sep="")])
                    } else {
                        cat(g_calls[id == input$choose_variance,paste0("ref ",reference,"   user ",user_sample,"   orig ",cons,"\n",exon_intron,"   @ ",gen_coord,"\nfwd mut ",mut_peak_base_fwd,"  \tpeak% ",round(mut_peak_pct_fwd,digits=1),"  \tS/N ",round(mut_s2n_abs_fwd,digits=1),sep="")])
                    }
                }, error = function(er){
                    if(grepl("NAs introduced",er)) cat("type an integer number")
                })
            } else cat("")
        } else cat("load .abi/.ab1 file")
    })

    update_chosen_variances <- observe({
        input$execute_btn
        isolate({
            if(loading_processed_files() != "not") {
                g_calls[id==as.numeric(input$choose_variance)]$user_sample <<- input$change_peak
                g_calls[id==as.numeric(input$choose_variance)]$set_by_user <<- TRUE
            }
        })
    })

#    output$chosenCheckboxes <- reactive({
#        return(loading_processed_files() != "not" & !is.null(g_choices))
#    })
#    outputOptions(output, "chosenCheckboxes", suspendWhenHidden = F)

#     output$table2 <- DT::renderDataTable({
#         input$execute_btn
#         input$delete_btn
#         if(loading_processed_files() != "not" & !is.null(g_choices)) {
#             DT::datatable(
#                 g_choices,
#                 rownames = checkboxRows(g_choices), escape = -1,
#                 options = list(dom = 'ti')
#             )
#         }
#     })

#     variace_selected_2 <- observe({
#         if(loading_processed_files() != "not") {
#             g_selected <<- g_choices[input$table2_selected,id]
#             #cat(input$table2_selected)
#         }
#     })

#     add_checkboxes <- reactive({
#         if(nrow(g_choices) > 0){
#             input$execute_btn
#             input$reset_btn
#             checkboxes <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '"',"")
#             for(i in 1:nrow(g_choices)) {
#                 if(g_choices[i]$id %in% g_selected) checkboxes[i] <- paste0(checkboxes[i]," checked>","")
#                 else checkboxes[i] <- paste0(checkboxes[i],">","")
#             }
#             return(checkboxes)
#         }
#     })

    add_checkboxes <- function(){
        checkboxes <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '"',"")
            for(i in 1:nrow(g_choices)) {
                if(g_choices[i]$id %in% g_selected) checkboxes[i] <- paste0(checkboxes[i]," checked>","")
                else checkboxes[i] <- paste0(checkboxes[i],">","")
            }
        return(checkboxes)
    }

    output$chosen_variances_table <- shiny::renderDataTable({
        if(varcall())
        if(!is.null(g_choices) && nrow(g_choices) > 0){
            input$execute_btn
            input$reset_btn
            if(loading_processed_files() != "not" & !is.null(g_choices)) {

                #add_checkbox_buttons <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '">',"")
                #add_edit_buttons <- paste0('<a class="go-edit" href="" data-id="', g_choices$id, '"><i class="fa fa-crosshairs"></i></a>')
                add_edit_buttons <- paste0('<input type="button" class="go-edit" value="edit" name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
                add_zoom_buttons <- paste0('<input type="button" class="go-zoom" value="zoom" name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
                add_checkbox_buttons <- add_checkboxes()

                cbind(Pick=add_checkbox_buttons, Edit=add_edit_buttons, Zoom=add_zoom_buttons, g_choices[,list(id=id,user_sample,call,reference)])
            }
        } #else { output$infobox <- renderPrint({ cat("no variances") }) }
    }, options = list(dom = "t",orderClasses=c(-1,-2,-3), paging=F, columnDefs=list(list(targets=c("_all"), searchable=F),list(targets=c(0,1,2), orderable=F, title="")))
    , escape=c(-1,-2,-3)
    , callback =
        "function(table) {
            table.on('change.dt', 'tr td input:checkbox', function() {
                setTimeout(function () {
                    Shiny.onInputChange('rows', $(this).add('tr td input:checkbox:checked').parent().siblings(':nth-child(4)').map(function() {
                        return $(this).text();
                    }).get())
                }, 10);
            });
            table.on('click', '.go-edit', function(e) {
                e.preventDefault();
                $el = $(this);
                var id_data = $el.data('id');
                Shiny.onInputChange('goEdit', {
                    id: id_data
                });
            });
            table.on('click', '.go-zoom', function(e) {
                e.preventDefault();
                $el = $(this);
                var id_data = $el.data('id');
                Shiny.onInputChange('goZoom', {
                    id: id_data
                });
            });
        }"
    )

    output$call_table <- shiny::renderDataTable({
    	if(loading_processed_files() != "not" & !is.null(g_calls)) { g_calls }
    }, options = list(paging=F, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    output$intens_table <- shiny::renderDataTable({
        if(loading_processed_files() != "not" & !is.null(g_intens)) { g_intens }
    }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    output$intens_table_rev <- shiny::renderDataTable({
        if(loading_processed_files() != "not" & !is.null(g_intens_rev)) { g_intens_rev }
    }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title=""))))

    variance_select <- observe({
        if(loading_processed_files() != "not") {
            g_selected <<- str_trim(input$rows)
            g_selected_zoom_index <<- 0
#            g_selected_edit_index <<- 0
        }
    })

    goEdit_handler <- observe({
        if(loading_processed_files() != "not") {
            if(is.null(input$goEdit)) return()
            updateTextInput(session,"choose_variance",value=paste0(input$goEdit$id))
        }
    })

    goZoom_handler <- observe({
        if(loading_processed_files() != "not") {
            if(is.null(input$goZoom)) return()
            #cat(paste0(input$goZoom$id,","))
            session$sendCustomMessage(type = 'zoom',message = paste0(g_calls[id==input$goZoom$id]$trace_peak))
        }
    })

    s2n_slider_handler <- observe({
        if(loading_processed_files()!= "not"){
            session$sendCustomMessage(type = 'mut_min',message = paste0(input$s2n_min))
        }
    })

    goClick_handler <- observe({
      if(loading_processed_files() != "not") {
        if(is.null(input$posClick)) return()
        updateTextInput(session,"choose_variance",value=paste0(input$posClick$id))
      }
    })

    reset_handler <- observe({
        input$reset_btn
        isolate({
            if(length(g_selected) != 0) {
                #g_choices <<- g_choices[-match(as.numeric(g_selected),gid)]
                for(i in as.numeric(g_selected)) {
                    g_calls[id==i,]$user_sample <<- g_calls[id==i,]$reference
                    g_calls[id==i,]$set_by_user <<- TRUE
                }
                #g_selected <<-  NULL
            }
        })
    })

    set_opacity <- observe({
        if(loading_processed_files() != "not"){
            session$sendCustomMessage(type = "opac_f",message = paste0(input$opacity_fwd/100))
            session$sendCustomMessage(type = "opac_r",message = paste0(input$opacity_rev/100))
        }
    })

})
