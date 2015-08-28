require(shiny)
library(sangerseqR)
source("helpers.R")

g_call      <- NULL             #annotated basecall data
g_ins       <- NULL             #intensities file
g_helperdat      <- NULL             #meta data used in graphs 
g_choices   <- NULL
g_selected  <- NULL
g_selected_zoom_index <- 0

#
g_abif <- NULL

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
            withProgress(message = paste('Loading file',file_name,'...',sep=" "), value = 1, {
                # TODO - make sure shiny only allows ab1 files (???)
                #        - otherwise chceck if file has correct format
                g_abif <- sangerseqR::read.abif(file)
                ins<- get_intensities(g_abif@data)
                res <-get_call_data(g_abif@data) 
                call <- res$call
                res$helperdat$max_y <-max(ins)
                res$helperdat$max_x <- nrow(ins)
                g_helperdat<<- res$helperdat
                call.dt <- data.table(call,key="id")
                g_call<<- call.dt        
                g_ins <<- ins
                
            })
            
            return("loaded")
        }
        return("not")
    })
    
    
    output$plot <- renderChromatography({
        if(loading_processed_files() != "not") {
#            withProgress(message="Rendering plot ...", value=1, {
                chromatography(g_ins,g_helperdat)
                #print(str(session.request))
#            })
        }
    })

    
    output$variance_info <- renderPrint({
        if(loading_processed_files() != "not") {
            if(input$choose_variance != "") {
                tryCatch({
                    cat(g_call[id == as.integer(input$choose_variance),paste("chosen variance: ",get("reference"),"(",get("exon_intron"),":",gen_coord,") -> ",call," with quality ",quality,sep="")])
                }, error = function(er){
                    if(grepl("NAs introduced",er)) cat("you may type in just integers")
    #                else cat("Some error")
                })
            } else cat("choose variance by number of peak (above)")
        } else cat("you must load some file first")
    })

    
    first_update_chosen_variances <- observe({
        if(loading_processed_files() != "not") {
          #TO DO: more sophisticated rules (might need to take intesitied into account)
          g_choices <<- g_call[call != reference  & seq_trim != "low_qual" & trace_peak != "NA"& !is.na(gen_coord)]         
        }
    })
   

    update_chosen_variances <- observe({
        input$execute_btn
        isolate({
            if(loading_processed_files() != "not") {
                if(is.null(g_choices))
                    g_choices <<- g_call[as.integer(input$choose_variance)]
                else {
                    if(!as.integer(input$choose_variance) %in% g_choices$id) {
                        new_variance <- g_call[as.integer(input$choose_variance)]
                        new_variance$call <- input$change_peak
                        g_choices <<- rbind(g_choices,new_variance)
                    } else if(g_choices[id == as.integer(input$choose_variance)]$call != input$change_peak)
                        g_choices[id == as.integer(input$choose_variance)]$call <<- input$change_peak
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
            add_edit_buttons <- paste0('<input type="button" class="go-edit" value="Edit" name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
            add_zoom_buttons <- paste0('<input type="button" class="go-zoom" value="Zoom" name="btn',g_choices$id,'" data-id="',g_choices$id,'"',">")
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
            cat(paste0(input$goZoom$id,","))
            session$sendCustomMessage(type = 'zoom_message',message = paste0(input$goZoom$id))
            
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