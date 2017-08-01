library(shiny)
library(data.table)
library(sangerseqR) #bioclite
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
    g_single_rev            <- NULL
    g_intrexdat             <- NULL             #intrex data used in graphs
    g_glassed_ref           <- NULL
    g_glassed_cod           <- NULL
    g_custom_ref            <- NULL
    g_custom_cod            <- NULL
    g_glassed_snp           <- NULL
    g_alignTo_options       <- c("TP53","ATM","NOTCH1","CALR")
    g_alignTo_description   <- list("TP53"   = "<b>TP53</b> description"
                                  , "ATM"    = "<b>ATM</b> description"
                                  , "NOTCH1" = "<b>NOTCH1</b> description"
                                  , "CALR"   = "<b>CALR</b> description")
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
    g_brush_fwd             <- 0
    g_brush_rev             <- 0
    g_not_loaded            <- ""
    g_reactval              <- reactiveValues()
    g_reactval$updateVar    <- 0
    g_refs_avail            <<- c("-","TP53","NOTCH1","ATM","CALR","Custom")
    #g_files = a variable that represents the data table containing individual samples and information about them
    TP53_demo               <- list(FWD_name=c("TP53 ; fwd ; low freq w frameshift"),
                                    FWD_file=c("data/demo/TP53/3low_freq_fsF.ab1"),
                                    REV_name=c("TP53 ; rev ; low freq w frameshift"),
                                    REV_file=c("data/demo/TP53/3low_freq_fsR.ab1"),
                                    REF=c("TP53"),
                                    mut_min=20,qual_thres_to_call=0,s2n_min=2,show_calls_checkbox=F,join_traces_checkbox=F,max_y_p=100,opacity=0,incorporate_checkbox=F,loaded=F,
                                    calls = "",
                                    status="new",
                                    id=1,
                                    brush_fwd_start = 0,
                                    brush_fwd_end = 0,
                                    brush_rev_start = 0,
                                    brush_rev_end = 0,
                                    coding = '',
                                    protein = '',
                                    VAF = '',
                                    dbSNP = '',
                                    dbSNP_id = '')
    g_files                 <- data.table(FWD_name=character(),FWD_file=character(),REV_name=character(),REV_file=character(),REF=character(),mut_min=numeric(),
                                           qual_thres_to_call=numeric(),s2n_min=numeric(),show_calls_checkbox=logical(),join_traces_checkbox=logical(),max_y_p=numeric(),opacity=numeric(),
                                           incorporate_checkbox=logical(),loaded=logical(),calls=character(),status=character(),id=integer(),brush_fwd_start=numeric(),
                                           brush_fwd_end=numeric(),brush_rev_start=numeric(),brush_rev_end=numeric(),coding=character(),protein=character(),VAF=character(),
                                           dbSNP=character(),dbSNP_id=character())
    #g_files[,`:=` ]

#     get_file <- reactive({
#        # if (!is.null(input$select_file)) return(input$select_file$datapath)
#         if (!is.null(input$select_file)) return(input$select_file)
#         else return(NULL)
#     })
    ex <- NULL
    btn_counter <- 0
    show_btn_counter <- 0

    #SAMPLE BROWSER STUFF

    output$samples_table <- DT::renderDataTable({

        updateRefs()
        loadSamples()
        goRef_handler()
        goDelete_sample_handler()
        goSwap_sample_handler()
        goChangeRev_handler()
        loading_processed_files()
        input$goLock
        loadDemo()



        disabled <- rep(FALSE,nrow(g_files))
        disabled[1] <- FALSE
        g_files[,id:= 1:nrow(g_files)]
        g_files <<- g_files
        if(nrow(g_files)==0){
            shinyjs::hide("export_btn")
            return(NULL)
        }
        add_load_buttons     <- shinyInput(actionButton, 1:nrow(g_files), 'loadSample_', label = NULL, onclick = 'Shiny.onInputChange(\"goLoadSamples\",  this.id + (Math.random()/10))',ico=rep("play",nrow(g_files)),class="btn btn-info" )
        add_delete_buttons   <- shinyInput(actionButton, 1:nrow(g_files), 'delSample_',  label = NULL, onclick = 'Shiny.onInputChange(\"goDeleteSamples\",  this.id + (Math.random()/10))',ico=rep("close",nrow(g_files)) ,dsbl = disabled,class="btn dlt_btn")
        add_swap_buttons     <- shinyInput(actionButton, 1:nrow(g_files), 'swapSample_', label = NULL, onclick = 'Shiny.onInputChange(\"goSwapSamples\",  this.id + (Math.random()/10))',ico=rep("exchange",nrow(g_files))  ,dsbl = disabled)
        add_reference_dropdown <- shinyInput(selectInput, 1:nrow(g_files),'selectGene_', choices=c("-",g_alignTo_options), selected = g_files[,REF],width="80px")
        add_reverse_dropdown <- shinyInputRev(selectInput,1:nrow(g_files),'chooseRev_',g_files,width="240px")
        out <- cbind(" "=add_delete_buttons,g_files[,list("forward"=FWD_name)],"swap"=add_swap_buttons,"<div title='Use dropdown menu to pair or unpair samples. Only unpaired reverse files appear here.' >reverse [?] </div>"=add_reverse_dropdown,"<div title='' >reference</div>"=add_reference_dropdown," "=add_load_buttons,g_files[,list("<div title='Confirmed (locked) variants appear here.'>status [?]</div>"=status,coding,protein,VAF,dbSNP)])
        table_out <- DT::datatable(out,escape=FALSE,selection = "none",
                                   style ="bootstrap",
                                   class="compact samplesdt",
                                   #id = "samplesdt",
                                   options=list("ordering"=FALSE,"paging"=FALSE,"searching"=FALSE,"autoWidth"=FALSE,"bInfo"=FALSE,"id"="samplesdt",
                                                initComplete = JS('function(settings, json) {
                                                                        $(\'[id*="selectGene"]\').change(function() {
                                                                                var sel = $(this).find(":selected").text();
                                                                                Shiny.onInputChange("goChangeRef",{id: this.id,
                                                                                                                   gene: sel});
                                                                        });
                                                                        $(\'[id*="chooseRev"]\').change(function() {
                                                                                var sel = $(this).find(":selected").text();
                                                                                Shiny.onInputChange("goChangeRev",{id: this.id,
                                                                                                                   name: sel});
                                                                                /*alert("changing rev "+this.id+" to "+ sel);*/
                                                                        });
                                                                  }')
                                                )
                                   )
        if(any(g_files$status!="new" && g_files$status!="error")){
            shinyjs::show("export_btn")
        }else{
            shinyjs::hide("export_btn")
        }
        return(table_out)
    })

    #Handlers for the Sample Browser

    loadSamples <- reactive({
        if(!is.null(input$browser_files)){
            g_not_loaded <- ""
            loaded <- ""
            error = ""
            ret <- samplesLoad(input$browser_files,output,g_files,g_alignTo_options,g_custom_ref)
            g_files <<- rbind(g_files[,!c("id"),with=FALSE],ret$loaded,fill=TRUE)
            g_files[,id:= 1:nrow(g_files)]
            g_files <<- g_files
            if(length(ret$not_loaded)>1){
                output$files <- renderPrint(paste0("<pre>Unable to load the following files: ",paste(ret$not_loaded[2:length(ret$not_loaded)],collapse="<br>"),"</pre>" ))
            }
        }
    })

    goChangeRev_handler <- reactive({
        if(!is.null(input$goChangeRev)){
            pos_at <- as.numeric(strsplit(input$goChangeRev$id,"_")[[1]][2])
            name <- as.character(input$goChangeRev$name)
            if(name=="-"){  #split
                rev_name <- g_files[pos_at]$REV_name
                rev_file <- g_files[pos_at]$REV_file
                ref      <- g_files[pos_at]$REF

                g_files[pos_at]$REV_name <- "-"
                g_files[pos_at]$REV_file <- "-"
                g_files[pos_at]$calls <- ""
                g_files[pos_at]$loaded <- FALSE
                g_files[pos_at]$status <- "new"
                g_files[pos_at]$brush_rev <- 0
                g_files[pos_at]$calls = ""
                g_files <<- rbind(g_files,c(list(FWD_name="-",FWD_file="-",
                                                 REV_name=rev_name,REV_file=rev_file,
                                                 REF=ref,id=nrow(g_files)+1),
                                                 mut_min=20,qual_thres_to_call=0,
                                                 s2n_min=2,show_calls_checkbox=F,
                                                 join_traces_checkbox=F,max_y_p=100,
                                                 opacity=0,incorporate_checkbox=F,
                                                 calls = "",loaded= FALSE,
                                                 status = "new",
                                                 brush_fwd_start = 0,brush_fwd_end=0,
                                                 brush_rev_start=0, brush_rev_end = 0,
                                                 coding = '', protein = '', VAF = '',dbSNP = '',dbSNP_id = ''))
            }else{          #combine
                setkey(g_files,id)
                rev_name <- name
                rm_id <- g_files[g_files[,REV_name == name]]$id
                rev_file <- g_files[g_files[,REV_name == name]]$REV_file
                pos_at <- as.numeric(strsplit(input$goChangeRev$id,"_")[[1]][2])
                g_files[pos_at]$REV_name <- rev_name
                g_files[pos_at]$REV_file <- rev_file
                g_files<<-g_files[!g_files[id==rm_id]]

            }
        }
    })
    goRef_handler <- reactive({
        if(!is.null(input$goChangeRef)){
            ref_id <- as.numeric(strsplit(input$goChangeRef$id,"_")[[1]][2])
            g_files[ref_id]$REF <-  as.character(input$goChangeRef$gene)
            g_files <<- g_files
        }
    })

    goDelete_sample_handler <- reactive({
        if(!is.null(input$goDeleteSamples)){
            isolate({
                delete_id <- floor(as.numeric(strsplit(input$goDeleteSamples, "_")[[1]][2]))/10
                #g_files <<- g_files[-delete_id]
                g_files <<- g_files[!g_files[,id==delete_id]]
                #js$delRow(delete_id)
            })
        }
    })

    goSwap_sample_handler <- reactive({
        input$goSwapSamples
        if(!is.null(input$goSwapSamples)){
            isolate({
                swap_id <- floor(as.numeric(strsplit(input$goSwapSamples, "_")[[1]][2]))/10
                swp_name <- g_files[id==swap_id]$FWD_name
                swp_file <- g_files[id==swap_id]$FWD_file
                g_files[id==swap_id]$FWD_name <- g_files[id==swap_id]$REV_name
                g_files[id==swap_id]$FWD_file <- g_files[id==swap_id]$REV_file
                g_files[id==swap_id]$REV_file <- swp_file
                g_files[id==swap_id]$REV_name <- swp_name
                g_files <<- g_files
                #js$swapRow(swap_id)

                updateSelectizeInput(
                    session, 'selectizeInput_1',server = FALSE,
                    options = list ("maxInput" = 5)
                )
            })
        }
    })

    #save settings in sample browser
    controls_listener <- observe({
        input$mut_min
        input$qual_thres_to_call
        input$max_y_p
        input$s2n_min
        input$show_calls_checkbox
        input$join_traces_checkbox
        input$max_y_p
        input$opacity
        input$incorporate_checkbox
        g_reactval$updateVar
        input$trim_fwd_end
        input$trim_fwd_start
        input$trim_rev_end
        input$trim_rev_start

        if(nrow(g_files[loaded==TRUE,]) ==1){
            g_files <<- g_files[loaded==TRUE,`:=`(mut_min=input$mut_min,
                                      qual_thres_to_call=input$qual_thres_to_call,
                                      s2n_min=input$s2n_min,
                                      show_calls_checkbox=input$show_calls_checkbox,
                                      join_traces_checkbox=input$join_traces_checkbox,
                                      max_y_p=input$max_y_p,
                                      opacity=input$opacity,
                                      incorporate_checkbox=input$incorporate_checkbox,
                                      brush_fwd_start=input$trim_fwd_start,
                                      brush_rev_end=input$trim_rev_end)]
        }

    })

    loading_processed_files <- reactive ({

        calls <- structure("error_reading_Rbin",class = "my_UI_exception")

        if(!is.null(input$goLoadSamples)){
            isolate({
                load_id <- floor(as.numeric(strsplit(input$goLoadSamples, "_")[[1]][2]))/10
                #save previously loaded
                if(nrow(g_files[loaded==TRUE,])==1){
                    tmpfile <- tempfile("calls")
                    save(g_calls,file=tmpfile)
                    g_files[loaded==TRUE,calls:=tmpfile]
                }
                g_files[,loaded:=F]
                g_files[load_id,loaded:=T]
                updateSliders(session,g_files)
                g_files<<-g_files
            })
        }

        single_rev <- FALSE
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
            if(ref=="Custom"){
                g_glassed_ref <<- g_custom_ref
                g_glassed_cod <<- g_custom_cod
            }else{
                g_glassed_ref <<- paste("data/refs/",ref,".glassed.intrex.fasta",sep="")
                g_glassed_cod <<- paste("data/refs/",ref,".glassed.codons.rdata",sep="")
            }

            snp_file      <-  paste("data/refs/",ref,".avsnp147.annovar.revcom.csv",sep="")
            if(file.exists(snp_file)){
                g_glassed_snp <<- fread(snp_file,header=FALSE)
            }else{
                g_glassed_snp <<- NULL
            }

            if (fwd_file_name == "-"){
                single_rev <- TRUE
                fwd_file <- rev_file
                rev_file <- NULL
            }
            else
                single_rev <- FALSE
            if (rev_file_name == "-")
                rev_file <- NULL

            g_single_rev <<- single_rev

            #if(substr(base,nchar(base),nchar(base))=="R"){
            if(single_rev){
                isolate({
                    files_info <<- paste0("fwd (F): ",fwd_file_name,"\nrev (R): ",rev_file_name," \n<em>aligned to: ",ref,"</em>",sep="")
                    #single_rev <- TRUE
                })
            }else if(is.null(ex)){
                isolate({
                    files_info <<- paste0("fwd (F): ",fwd_file_name,"\nrev (R): ",rev_file_name," \n<em>aligned to: ",ref,"</em>",sep="")
                })
            }

            withProgress(message = paste('loading abi file...',sep=" "), value = 1, {

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
                    error = function(e){
                        output$files <- renderPrint(paste0("<pre>error while loading calls from abi file : ",e$message,"</pre>" ))
                        g_files[loaded==T,status:="<font color='red'>error</font>"]
                    }
                )

                if(!is.null(called)){
                    g_minor_het_insertions  <<- NULL
                    g_stored_het_indels     <<- list()
                    g_indels_present        <<- FALSE
                    g_qual_present          <<- called$qual_present
                    intensified  <-  get_intensities(g_abif,g_abif_rev,calls=called$calls,deletions=called$deletions,norm=FALSE,single_rev)
                    g_intens     <<- intensified$intens
                    g_intens_rev <<- intensified$intens_rev
                    tryCatch(
                        calls        <-  annotate_calls(calls=intensified$calls,intens=intensified$intens,intens_rev=intensified$intens_rev,g_glassed_cod),
                        error = function(e){output$files <- renderPrint(paste0("<pre>error while loading calls from abi file : ",e$message,"</pre>" ))}
                    )
                    calls        <-  adjust_ref_mut(calls,g_intens_rev)
                    g_max_y      <<- max(c(max(g_intens[,list(A,C,G,T)]),if(is.null(g_intens_rev)) 0 else max(g_intens_rev[,list(A,C,G,T)])))
                    #intrex contains intesities coordinates of start and end of introns/exons with the sequence id (position in sequence coordinates)
                        intrexdat            <- list()
                        intrexdat$intrex     <- list()
                        intrexdat$intrex     <- setnames(calls[!is.na(exon_intron),list(max(id)-min(id)+1,min(trace_peak),max(trace_peak)),by = exon_intron],c("attr","length","trace_peak","end"))
                        intrexdat$intrex     <- setnames(merge(intrexdat$intrex,calls[,list(id,trace_peak)],by="trace_peak"),"trace_peak","start")
                        cs<- unlist(lapply(intrexdat$intrex$id,function(x){calls[id==x]$coding_seq}))
                        intrexdat$intrex[,coding_seq:=cs]
                        intrexdat$max_x      <- max(c(nrow(g_intens),nrow(g_intens_rev))) # these numbers should be the same
                        intrexdat$new_sample <- TRUE
                    g_intrexdat       <<- splice_variants(intrexdat)
                    calls             <-  data.table(calls,key="id")

                    # temporarily switching off functionality 2 sept 16, Karol
                    # g_noisy_neighbors <<- get_noisy_neighbors(calls)
                    if(!called$qual_present){
                        files_info <- paste0(files_info,HTML("\n<strong style=\"color: red;\">no Phred qualities!</strong>"))
                    }
                    files_info <- paste0("<pre>",files_info,"</pre>")
                    output$files      <-  renderPrint({cat(files_info)})
                    g_new_sample      <<- TRUE

                    #initialize or load trim brushes
                    lapply(c("trim_fwd_start","trim_fwd_end","trim_rev_start","trim_rev_end"),function(x){updateNumericInput(session,max = nrow(calls),inputId=x)})

                    if(g_files[loaded==TRUE,]$status =="new"){
                        g_files[loaded==TRUE,]$brush_fwd_start<<-calls[call!="-",][25]$id
                        g_files[loaded==TRUE,]$brush_fwd_end<<-calls[nrow(calls)-25]$id
                    }
                    updateNumericInput(session,value=g_files[loaded==TRUE,]$brush_fwd_start,inputId = "trim_fwd_start" )
                    updateNumericInput(session,value=g_files[loaded==TRUE,]$brush_fwd_end,inputId = "trim_fwd_end" )

                    if(!is.null(g_abif_rev)){
                        if(g_files[loaded==TRUE,]$status =="new"){
                            g_files[loaded==TRUE,]$brush_rev_start<<-calls[call_rev!="-",][25]$id
                            g_files[loaded==TRUE,]$brush_rev_end<<-calls[nrow(calls[call_rev!="-",])-25]$id
                        }
                        updateNumericInput(session,value=g_files[loaded==TRUE,]$brush_rev_start,inputId = "trim_rev_start" )
                        updateNumericInput(session,value=g_files[loaded==TRUE,]$brush_rev_end,inputId = "trim_rev_end" )
                    }else{
                        updateNumericInput(session,value=0,inputId = "trim_rev_start" )
                        updateNumericInput(session,value=nrow(calls),inputId = "trim_rev_end" )
                    }
                    #if(g_single_rev){
                    #    g_brush_fwd <<- calls[nrow(calls[call!="-",])-25]$trace_peak

                    #}
                    #if reloading a previously loaded file
                    if(nrow(g_files[loaded==T,]) == 1){                  #g_calls saved from previous session we test if they are compatible to reload
                        if(g_files[loaded==T,]$calls != ""){
                            load(g_files[loaded==T]$calls)
                        }
                    }

                    updateTabsetPanel(session,'tabs',selected = "main")

                } else {
                    return(structure("error_reading_Rbin",class = "my_UI_exception"))
                }
            })
        }
        g_stored_het_indels <<- list()
        g_calls <<- NULL
        return(calls)
    })

    get_mut_min <- eventReactive(input$mut_min,{
        return(input$mut_min)
    })

    varcall <- reactive({
        if(class(loading_processed_files())[1] != "my_UI_exception") {

            if(g_files[loaded==T,]$calls == "") g_files[loaded==T,status:="viewed"]
            updateChosenVariants()
            goReset_handler()
            g_reactval$updateVar
            #goBrush_fw()

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
            foo <- get_mut_min()
            g_calls <<- call_variants(g_calls,input$qual_thres_to_call,foo,input$s2n_min,g_stored_het_indels,g_brush_fwd,g_brush_rev,input$incorporate_checkbox,g_single_rev)
            #g_calls <<- call_variants(g_calls,input$qual_thres_to_call,foo,input$s2n_min,g_stored_het_indels,0,0,input$incorporate_checkbox,g_single_rev)
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
            g_choices <<- get_choices(g_calls,g_glassed_ref)
            return(TRUE)

        } else return(FALSE)
    })

    #
    # Render functions reacting to varcall
    #
    output$plot <- renderChromatography({
        if(varcall()) {
            g_intrexdat$max_y <- (g_max_y*100)/input$max_y_p
            withProgress(message = paste('Plotting chromatogram.',sep=" "), value = 1 ,{
                ret<-chromatography(g_intens,g_intens_rev,g_single_rev,g_intrexdat,
                                    g_calls,g_choices,g_new_sample,
                                    g_noisy_neighbors,
                                    input$show_calls_checkbox,input$show_qual_checkbox,
                                    g_qual_present,
                                    g_calls[input$trim_fwd_start]$trace_peak,g_calls[input$trim_fwd_end]$trace_peak,
                                    g_calls[input$trim_rev_start]$trace_peak,g_calls[input$trim_rev_end]$trace_peak)
                g_new_sample <<- FALSE
                return(ret)
            })
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
    updateChosenVariants <- reactive({
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
            #input$reset_btn
            #input$lo
            if(varcall() & !is.null(g_choices)) {
                g_view <<- getView(g_calls,g_choices,g_glassed_snp)
                g_view <<- applyFilters(g_view,input$trim_fwd_start,input$trim_fwd_end,input$trim_rev_start,input$trim_rev_end)
                #add_goto_buttons     <- shinyInput(actionButton, g_view$id, 'button_', label = "go", onclick = 'Shiny.onInputChange(\"goGoto\",  this.id+ (Math.random()/10))' )
                tp <- g_view$trace_peak
                add_goto_buttons     <- shinyInput(actionButton, g_view$id, 'button_', trace_peak = tp, label = "go", onclick = 'posClick(parseInt(this.id.split("_")[1]),parseInt(this.id.split("_")[2]));')
                add_reset_buttons    <- shinyInput(actionButton, g_view$id, 'button_', label = NULL, ico=rep("close",nrow(g_view)),onclick = 'Shiny.onInputChange(\"goReset\",  this.id)',class="btn dlt_btn" )
                add_lock_buttons     <- shinyInput(actionButton, g_view$id, 'button_', label = NULL, onclick = 'Shiny.onInputChange(\"goLock\", this.id+ (Math.random()/10));if($(this).children(":first").attr("class")=="fa fa-unlock"){$(this).children().addClass(\'fa-lock\').removeClass(\'fa-unlock\');}else{$(this).children().addClass(\'fa-unlock\').removeClass(\'fa-lock\');}',ico = unlist(lapply(g_view$set_by_user, function(x){if(isTRUE(x)){"lock"}else{ "unlock"}})),class="btn btn-success" )
                if(nrow(g_view)>0){
                    out<-cbind(" "=add_goto_buttons, " "=add_reset_buttons, "<div title='Confirmed variants are kept for the session even if you change parameters,\nand appear in the samples panel from where they can be exported with the green export button.'>confirm [?]</div>"=add_lock_buttons, g_view[,list("call position (start)"=id,"genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein,"pri peak %"=sample_peak_pct,"sec peak %"=mut_peak_pct)])
                    if(!is.null(g_glassed_snp)){
                       out <- cbind(out,g_view[,list("dbSNP"=dbSNP)])
                    }
                }else{
                    out <- g_view[,list("call position (start)"=id,"genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein,"pri peak %"=sample_peak_pct,"sec peak %"=mut_peak_pct)]
                }
                tableout<-DT::datatable(out
                                        , escape=FALSE
                                        #, class = "compact"
                                        , selection = "none"
                                        , style= 'bootstrap'
                                        , options=list("paging"=FALSE,"searching"=FALSE,"ordering"=FALSE,"autoWidth"=FALSE)
                                        )
            }
        }
    }
    #, options = list(dom = "t",orderClasses=c(-1,-2,-3,-4), paging=F, columnDefs=list(list(targets=c("_all"), searchable=F),list(targets=c(0,1,2,3), orderable=F, title="")))

    )

#    variant_select <- observe({
#        if(varcall()) {
#            g_selected <<- str_trim(input$rows)
#            g_selected_goto_index <<- 0
#        }
#    })

    #
    # data table buttons  ### replaced with a direct call to chromatography.js
    #
    #goGoto_handler <- observe({
    #    if(is.null(input$goGoto)) return()
    #    goto_id <- floor(as.numeric(strsplit(input$goGoto, "_")[[1]][2]))/10
        #session$sendCustomMessage(type = 'goto',message = paste0(g_calls[id==goto_id]$trace_peak))
    #    updateTextInput(session,"choose_call_pos",value=paste0(goto_id))
    #})

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
            lock_id <- strsplit(input$goLock, "_")[[1]][2]

            if(length(strsplit(lock_id,".",fixed=TRUE)[[1]]) >2){
                lock_id <- strsplit(lock_id,".",fixed=TRUE)[[1]]
                lock_id <- as.numeric(paste0(lock_id[1],".", lock_id[2]))
            }else{
                lock_id <- floor(as.numeric(lock_id))/10
            }


            g_view[id==lock_id]$set_by_user <<- !g_view[id==lock_id]$set_by_user
            coding  <- g_view[id==lock_id]$coding
            updateTextInput(session,"choose_call_pos",value=paste0(lock_id))
                    if(length(grep("ins|del|dup",coding)) > 0){ #locking indels
                        len <- length(unique(na.omit(as.numeric(unlist(strsplit(unlist(coding), "[^0-9]+"))))))
                        if(is.null(g_stored_het_indels[[coding]])){
                            g_stored_het_indels[[coding]] <<- lock_id
                            for(i in 0:(len-1)){
                                g_calls[id==(lock_id + i)]$set_by_user <<- !g_calls[id==(lock_id +i)]$set_by_user
                            }
                        }else{
                            g_stored_het_indels[[coding]] <<- NULL
                            for(i in 0:(len-1)){
                                g_calls[id==(lock_id+i)]$set_by_user <<- !g_calls[id==(lock_id +i)]$set_by_user
                            }
                        }

                    } else {     #locking SNPs
                g_calls[id==lock_id]$set_by_user <<- !g_calls[id==lock_id]$set_by_user
            }

            if(nrow(g_view[set_by_user==TRUE])>0){
                #prots <- paste(g_view[set_by_user == TRUE]$protein,collapse="")
                #g_files<<-g_files[loaded==TRUE,status:=paste0("<b>confirmed</b>: ",paste(g_view[set_by_user == TRUE]$coding,collapse=";"),prots)]
                n = nrow(g_view[set_by_user==TRUE])

                g_files<<-g_files[loaded==TRUE,':='(status   = paste0(n," mut. found",collapse=""),
                                                    coding   = paste(g_view[set_by_user == TRUE]$coding,  collapse="<br>"),
                                                    protein  = paste(g_view[set_by_user == TRUE]$protein, collapse="<br>"),
                                                    dbSNP    = paste(g_view[set_by_user == TRUE]$dbSNP,   collapse="<br>"),
                                                    dbSNP_id = paste(gsub("<a .*>(.*)</a>","\\1",g_view[set_by_user == TRUE]$dbSNP) ,collapse="<br>"),
                                                    VAF      = paste(g_view[set_by_user == TRUE]$VAF,     collapse=" <br> ")
                                                    )]
            }else{
                g_files<<-g_files[loaded==TRUE,status:="viewed"]
            }
            #output$goLock <- renderUI({actionButton(input$goLock, icon = icon("lock"))})
        })
        return(T)
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

    goBrush_fw <- observe({
        if(is.null(input$brush_fw)) return()
        last <- g_brush_fwd
        g_brush_fwd <<- input$brush_fw$coord
        #print(input$brush_fw$coord)
        #print(g_brush_fwd)
        if(last <= g_brush_fwd){
            if(nrow(g_choices[trace_peak > last][trace_peak < g_brush_fwd])==0)
                return()
        }
        g_reactval$updateVar <- runif(1,0,1)
    })
    goBrush_rv <- observe({
        if(is.null(input$brush_rv)) return()
        last <- g_brush_rev
        g_brush_rev <<- input$brush_rv$coord
        if(g_single_rev){     #we use brush_fwd (as we only have forward calls) but visualised and logic as brush rev
            g_brush_rev <- last
            g_brush_fwd <<- input$brush_rv$coord
        }else{
            if(last>=g_brush_rev){
                if(nrow(g_choices[trace_peak_rev < last][trace_peak_rev>g_brush_rev])==0)
                   return()
            }
        }
        g_reactval$updateVar <- runif(1,0,1)
    })
    #output$helpButton <- renderUI({
    #    actionButton("toggle_help",icon = icon("question"), label= label())
    #})
    observe({
        if(!is.null(input$toggle_help )){
            #show_btn_counter <<- show_btn_counter + 1
            if(input$toggle_help %% 2 == 0) html("toggle_help",'<i class="fa fa-question"></i> hide help')
            else html("toggle_help",'<i class="fa fa-question"></i> show help')
        }
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

    #######################
    # alignTo (start) #
    #######################

    #g_alignTo_options <- c("TP53","ATM","NOTCH1","CALR")
    g_alignTo_description <- list("TP53"   = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/NM_000546.5' target='_blank'>NM_000546.5</a> <br> GRCh38"
                                , "ATM"    = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/NM_000051.3' target='_blank'>NM_000051.3</a> <br> hg19"
                                , "NOTCH1" = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/NM_017617.4' target='_blank'>NM_017617.4</a> <br> GRCh38"
                                , "CALR"   = "<a href='https://www.ncbi.nlm.nih.gov/nuccore/NM_004343.3' target='_blank'>NM_004343.3</a> <br> GRCh38")

    output$alignTo_new <- renderUI({
        updateRefs()
        setCustomRef()
        UI_out <- lapply(g_alignTo_options, function(x) paste0(div(id = paste0(x, '_alignTo_div')
                                                                 , class = 'alignTo_class'
                                                                 , style = "margin: 0.3em;
                                                                            padding: 0.4em;
                                                                            height: 8.0em;
                                                                            min-width: 90px;
                                                                            background-color: #F5F5F5;
                                                                            box-shadow: 0.1em 0.1em 0.3em #888888;"
                                                                 , h4(x)
                                                                 , div(HTML(g_alignTo_description[[x]]))
                                                                 , if(x=="TP53") actionLink(paste0(x, "_alignTo_lnk"), label = "[load demo]") else "-"
                                                                 )
                                                               )
                         )
        UI_out <- paste0(UI_out, collapse = "</td><td>")
        UI_out <- paste0("<table><tr><td>", UI_out, "</td></tr></table>")
        return(HTML(UI_out))
    })


    loadDemo <- reactive({
        if(!is.null(g_reactval$link)){
            if(g_reactval$link=="TP53"){
                g_files <<- rbind(g_files,TP53_demo)
            }
        }
    })

    updateRefs <-function() {
        #loadGBK()
        #setCustomRef()
        g_alignTo_options <<- input$additionalRefs
        if(!is.null(g_custom_ref)){
            g_alignTo_options <<- c(g_alignTo_options,"Custom")
        }

        lapply(g_alignTo_options, function(x) {
            observeEvent(input[[paste0(x, "_alignTo_lnk")]], {
                g_reactval$link <- x
            })
        })
    }

    loadGBK <- observeEvent(input$custom_gb,{
        #if(!is.null(input$custom_gb)){
            withProgress(
                tryCatch({
                        res <- get_gbk_info(session,input$custom_gb)
                        choices<-NULL
                        for(i in 1:length(res$CDS)){
                            option <- paste0("GENE: ",res$CDS[[i]]$gene,", PRODUCT: ",res$CDS[[i]]$product," (",nchar(res$CDS[[i]]$translation),"aa)",collapse = "")
                            choices[[option]] <- i
                        }

                        showModal({modalDialog(
                                #p('Found ',length(res$mRNA),' "mRNA" features in GBK file '),
                                p('searching for "CDS" features in the GenBank file...'),
                                #if (failed)
                                #div(tags$b("Invalid name of data object", style = "color: red;")),

                                radioButtons("radio_gene_sel", label = h3("choose one of the following gene products"),
                                             choices = choices,
                                             selected = 1,width = "100%"),
                                footer = tagList(
                                    modalButton("cancel"),
                                    actionButton("gbk_ok", "OK")
                                )
                        ,size='m')})
                        #ret <- process_gbk(session,input$custom_gb)
                        #g_custom_cod <<- ret$custom_cod
                        #g_custom_ref <<- ret$custom_fasta
                    }
                    ,error = function(e){
                        output$files <- renderPrint(paste0("<pre>error while reading GenBank file: ",e$message,"</pre>" ))
                        g_custom_cod <<- NULL
                        g_custom_ref <<- NULL
                    }
                ),message = "processing GenBank file..."
            )
        #}
    })

    setCustomRef <- reactive({
        #if(!is.null(input$custom_gb) && !is.null(input$gbk_ok) && !(input$gbk_ok==0)){
        if( !is.null(input$gbk_ok) &&!(input$gbk_ok==0)){
        #input$gbk_ok
        isolate({
            if(!is.null(input$custom_gb)){
                withProgress(
                    tryCatch({
                        removeModal()
                        showModal({modalDialog(
                            #p('found ',length(res$mRNA),' "mRNA" features in GBK file '),
                            p('extracting sequence information from GenBank file...'),
                            #if (failed)
                            #div(tags$b("invalid name of data object", style = "color: red;")),

                            footer = tagList(
                            )
                            ,size='m')})
                        ret <- process_gbk(session,input$custom_gb,as.numeric(input$radio_gene_sel))
                        g_custom_cod <<- ret$custom_cod
                        g_custom_ref <<- ret$custom_fasta
                        removeModal()
                        g_alignTo_options <<- input$additionalRefs
                        if(!is.null(g_custom_ref)){
                            g_alignTo_options <<- c(g_alignTo_options,"Custom")
                        }

                        lapply(g_alignTo_options, function(x) {
                            observeEvent(input[[paste0(x, "_alignTo_lnk")]], {
                                g_reactval$link <- x
                            })
                        })

                    }
                    ,error = function(e){
                        removeModal()
                        output$files <- renderPrint(paste0("<pre>error while reading GenBank file: ",e$message,"</pre>" ))
                        g_custom_cod <<- NULL
                        g_custom_ref <<- NULL
                    }
                    ),message = NULL,detail=NULL,style="old"
                )
            }
        })
        }
    })

    #####################
    # alignTo (end) #
    #####################


    #EXPORT

    output$export_btn <- downloadHandler(
        filename = function() {
            paste('data-', Sys.Date(), '.xlsx', sep='')
        },
        content = function(con) {
            out<-data.table("sample"=character(),"coding variant"=character(),"protein variant"=character(),"VAF"=character(),"dbSNP"=character())

            for(i in 1:nrow(g_files)){
                if(g_files[i]$status=="new"){
                    var = "NA"
                }else if(g_files[i]$status == "viewed"){
                    var = "wt"
                }else if(g_files[i]$status == "error"){
                    next
                }
                else{
                    var = g_files[i]$coding
                }
                out <- rbind(out,list(paste0(g_files[i]$FWD_name,":",g_files[i]$REV_name),
                                      gsub("<br>",";", var),
                                      gsub("<br>",";", g_files[i]$protein),
                                      gsub("<br>",";", g_files[i]$VAF),
                                      gsub("<br>",";", g_files[i]$dbSNP_id) ))
            }
                #out<-rbind(out,list("1","2","2"))


            #g_selected <- g_view$id[as.numeric(input$chosen_variants_table_rows_selected)]

            #for(i in 1:nrow(g_view)) {
            #    if(g_view[i]$id %in% g_selected) out<-rbind(out,g_view[i,list("genomic coordinate"=gen_coord,"coding variant"=coding,"protein variant"=protein)])
            #}
            write.xlsx(out,con)
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
    output$upload_file <- renderUI({
        input$img_help
        images <- data.table(
                   c(
                      "select ref from list"
                     ,"load custom GenBank ref"
                     ,"file upload"
                     ,"scroll"
                     ,"zoom"
                     ,"export detected variants"
                     )
                  ,c(
                      'src = "select_ref_from_list.gif" style="width:640px;height:394px;"'
                     ,'src = "custom_ref.gif" style="width:640px;height:394px;"'
                     ,'src = "LoadFile.gif" style="width:640px;height:394px;"'
                     ,'src = "navbar_scroll.gif" style="width:640px;height:394px;"'
                     ,'src = "navbar_zoom.gif"'
                     ,'src = "export.gif" style="width:960px;height:561px;"'
                     )
                  )

        UI_out <- paste0('<img ',images[V1==input$img_help]$V2,' ></img>')
        return(HTML(UI_out))
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
