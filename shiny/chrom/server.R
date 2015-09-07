library(shiny)
library(sangerseqR)
source("helpers.R")

g_calls      <- NULL             #annotated basecall data
g_intens    <- NULL             #intensities file
g_helperdat <- NULL             #helper data used in graphs
g_choices   <- NULL
makeReactiveBinding("g_choices")  #need someones opinion on this
g_selected  <- NULL
g_selected_zoom_index <- 0
g_chrom <- NULL
g_max_y <- NULL
g_abif  <- NULL

shinyServer(function(input,output,session) {

    get_file <- reactive({
       # if (!is.null(input$select_file)) return(input$select_file$datapath)
        if (!is.null(input$select_file)) return(input$select_file)
        else return(NULL)
    })

    loading_processed_files <- reactive ({
        
        if(!is.null(get_file())) {
            file <- get_file()$datapath
            file_name <- get_file()$name
            ret <- "not"
            withProgress(message = paste('...',sep=" "), value = 1, {
                # TODO - make sure shiny only allows ab1 files (???)
                #        - otherwise chceck if file has correct format
                g_abif      <- sangerseqR::read.abif(file)
                res <- ""
                tryCatch(res <- suppressWarnings(get_call_data(g_abif@data,input$rm7qual_thres,input$qual_thres,input$aln_min)),
                                                 error = function(e){output$variance_info <- renderPrint(paste0("An error occured while loading calls from abi file with the following error message: ", 
                                                                                                                e$message ))})
                if(res!=""){
                    output$variance_info <- renderPrint("")
                    intens      <- get_intensities(g_abif@data,norm=FALSE)     
                    calls        <- res$calls
                    g_max_y     <<- max(intens)
                    res$helperdat$max_x <- nrow(intens)
                    g_helperdat <<- res$helperdat
                    g_calls      <<- data.table(calls,key="id")
                    g_choices   <<- g_calls[user_mod != reference & user_mod != "low qual" & trace_peak != "NA" & !is.na(gen_coord)]
                    g_intens    <<- intens
                    ret<-"loaded"
                }
            })
            return(ret)
        }
        return("not")
    })
    #don't see the point of this anymore moved to loading processed files
#    first_update_chosen_variances <- observe({
#      if(loading_processed_files() != "not") {
#        #TO DO: more sophisticated rules (might need to take intensities into account)
#        g_choices <<- g_calls[user_mod != reference & user_mod != "low_qual" & trace_peak != "NA" & !is.na(gen_coord)]
#      }
#    })

    output$plot <- renderChromatography({
        if(loading_processed_files() != "not") {
            g_helperdat$max_y =  (g_max_y*100)/input$max_y_p
            chromatography(g_intens,g_helperdat,g_calls,g_choices,input$intens_guide_line)
        }
    })
    output$variance_info <- renderPrint({
        print("anything")
        input$pos_click
    })
    output$variance_info <- renderPrint({
        if(loading_processed_files() != "not") {
            if(input$choose_variance != "") {
                tryCatch({
                    cat(g_calls[id == input$choose_variance,paste("ref:",reference," call:",call," user:",user_mod," on ",exon_intron," at ",gen_coord," with quality ",quality,sep="")])
                }, error = function(er){
                    if(grepl("NAs introduced",er)) cat("type an integer number")
    #                else cat("Some error")
                })
            } else cat("")
        } else cat("load .abi/.ab1 file")
    })
    

    update_chosen_variances <- observe({
        input$execute_btn
        isolate({
            if(loading_processed_files() != "not") {
                if(is.null(g_choices))
                    g_choices <<- g_calls[input$choose_variance]
                else {
                    if(!input$choose_variance %in% g_choices$id) {
                        #new_variance <- g_calls[as.integer(input$choose_variance)]
                        new_variance <- g_calls[id==input$choose_variance]
                        new_variance$user_mod <- input$change_peak
                        g_choices <<- rbind(g_choices,new_variance)
                    } else if(g_choices[id == input$choose_variance]$user_mod != input$change_peak)
                        g_choices[id == input$choose_variance]$user_mod <<- input$change_peak
                }
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
#
#
#     variace_selected_2 <- observe({
#         if(loading_processed_files() != "not") {
#             g_selected <<- g_choices[input$table2_selected,id]
#             #cat(input$table2_selected)
#         }
#     })

    add_checkboxes <- reactive({
        input$execute_btn
        input$delete_btn
        checkboxes <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '"',"")
        for(i in 1:nrow(g_choices)) {
            if(g_choices[i]$id %in% g_selected) checkboxes[i] <- paste0(checkboxes[i]," checked>","")
            else checkboxes[i] <- paste0(checkboxes[i],">","")
        }
        return(checkboxes)
    })

    output$chosen_variances_table <- shiny::renderDataTable({
        input$execute_btn
        input$delete_btn
        if(loading_processed_files() != "not" & !is.null(g_choices)) {
#            add_checkbox_buttons <- paste0('<input type="checkbox" name="row', g_choices$id, '" value="', g_choices$id, '">',"")
            #add_edit_buttons <- paste0('<a class="go-edit" href="" data-id="', g_choices$id, '"><i class="fa fa-crosshairs"></i></a>')
            add_edit_buttons <- paste0('<input type="button" class="go-edit" value="edit" name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
            add_zoom_buttons <- paste0('<input type="button" class="go-zoom" value="zoom" name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
            add_checkbox_buttons <- add_checkboxes()
            cbind(Pick=add_checkbox_buttons, Edit=add_edit_buttons, Zoom=add_zoom_buttons, g_choices)
        }
    }, options = list(orderClasses=c(-1,-2,-3), paging=F, columnDefs=list(list(targets=c(0,1,2), searchable=F, orderable=F, title="")))
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
    	if(loading_processed_files() != "not" & !is.null(g_calls)) {
    		g_calls
    	}
    }, options = list(paging=F, columnDefs=list(list(searchable=F, orderable=F, title="")))
    )

    output$intens_table <- shiny::renderDataTable({
        if(loading_processed_files() != "not" & !is.null(g_intens)) {
            g_intens
        }
    }, options = list(paging=T, columnDefs=list(list(searchable=F, orderable=F, title="")))
    )

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

#     edit_handler <- observe({
#         input$edit_btn
#         isolate({
#             if(length(g_selected) != 0) {
#                 index <- g_selected_zoom_index + 1
#                 g_selected_zoom_index <<- (1 + g_selected_zoom_index) %% length(g_selected)
#                 updateTextInput(session,"choose_variance",value=paste0(g_selected[index]))
#             }
#         })
#     })

#     zoom_handler <- observe({
#         input$zoom_btn
#         isolate({
#             if(length(g_selected) != 0) {
#                 index <- g_selected_zoom_index + 1
#                 g_selected_zoom_index <<- (1 + g_selected_zoom_index) %% length(g_selected)
#                 cat(paste0(g_selected[index],","))
#                 session$sendCustomMessage(type = 'zoom_message',message = paste0(g_selected[index]))
#             }
#         })
#     })

    delete_handler <- observe({
        input$delete_btn
        isolate({
            if(length(g_selected) != 0) {
                g_choices <<- g_choices[-match(g_selected,id)]
                g_selected <<-  NULL
            }
        })
    })
})
