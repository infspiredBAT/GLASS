library(shiny)
library(shinyjs)
library(data.table)
library(rjson)
library(htmlwidgets)
library(chromatography)
library(stringr)
library(stringi)
library(DT)
library(shinyBS)
#source("JS.R")

shinyUI(
    fluidPage(
        useShinyjs(),
        theme = "simplex3.3.6.css", # http://bootswatch.com/ | sandstone/simplex/flatly/darkly
		tags$head(
		    includeCSS("www/samples.css"),
		    tags$head(HTML("<link href='https://fonts.googleapis.com/css?family=Inconsolata' rel='stylesheet' type='text/css'>")),
            tags$title("GLASS"),
		    tags$script('window.onbeforeunload = function() {
		                    return "You will need to start over if you navigate away!";
		                };'),
            tags$script('Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {
              		        var el = $("#" + x);
              		        el.replaceWith(el = el.clone(true));
              		        var id = "#" + x + "_progress";
              		        $(id).css("visibility", "hidden");
              		    });'),
            tags$style(HTML(".DTFC_LeftBodyLiner { width: 100% !important; }"))
		),

		fluidRow(
			column(1, HTML("&nbsp&nbsp<b><font size=6em>GLASS</font></b>")),
			column(1,HTML('<div  style="padding:.3em .5em"><a href="javascript:void(0)" id="app-disclaimer-link" onclick="$(\'#disclaimer-modal\').modal(\'show\')" >
                                ver 0.2.18 <br> (2017-Oct-5) 
                           </a></div>'), 
			       includeHTML("www/log.html")
			       ),
			column(10, HTML("<br><font size=2em>&nbsp&nbsp | R&D by <a href=http://bat.infspire.org target=_blank>bat.infspire.org</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.ceitec.eu/ceitec-mu/medical-genomics/rg34 target=_blank>Medical Genomics @ CEITEC MU</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.ericll.org target=_blank>European Research Initiative on CLL / ERIC</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.igcll.org target=_blank>IgCLL group</a></font>&nbsp&nbsp | IT by <a href=https://metavo.metacentrum.cz/en/ target='_blank'>CESNET/MetaCentrum</a>&nbsp&nbsp | @ <a href=mailto:bat@infspire.org target='_blank'>email us</a>"))
		),
		fluidRow(
			column(2, HTML(paste("&nbsp<i>assisted and standardised assessment</br>&nbsp&nbsp&nbspof gene variations from Sanger data</i></br></br>"))),
			column(9, htmlOutput("files"))
		),
		tabsetPanel(id = 'tabs',
            tabPanel('samples',value = 'smpl_brws',icon = icon("flask"),
				fluidRow(
					column(12,wellPanel(tags$div(HTML(paste("<div style=\"font-family:'Inconsolata';font-size:1.1em;\">
                                            <b>quick</b> &nbsp(0) hovering over '[?]'s will provide quick help tips </br>
                                            <b>guide</b> &nbsp(1) select / upload relevant references beforehand, GLASS will align and auto-orientate against them </br>&nbsp&nbsp&nbsp&nbsp&nbsp
                                                         &nbsp(2) upload ABI files from supported references, forward and/or reverse, and with unique names </br>&nbsp&nbsp&nbsp&nbsp&nbsp
                                                         &nbsp(3) after processing and auto-detection of most file properties, pair / unpair / swap / delete / change reference as necessary </br>&nbsp&nbsp&nbsp&nbsp&nbsp
                                                         &nbsp(4) click the blue 'play' button to load the file(s) </br>&nbsp&nbsp&nbsp&nbsp&nbsp
                                                         &nbsp(5) when done in 'variants' panel (find separate instructions there), confirmed variants will appear under 'status' and are exportable with the green 'export variants' button
                                        </div>"), sep = ""))))
					),
                fluidRow(
                    column(12,
                        tags$div(title="",HTML(paste("<hr><div style='display:inline;'><h4 style='display:inline;'>references</h4> &nbsp <b style='display:inline;' id ='ref_help'> [?]</b></div>")))
                        
                        # checkboxInput("showWP3", HTML("Manage References <i class='fa fa-wrench' aria-hidden='true'></i>")),
                        # conditionalPanel(condition="input.showWP3",
                             ,fluidRow(
                                 column(3,
                                        selectizeInput("additionalRefs", "select from our list of curated references",c("TP53","ATM","NOTCH1","CALR"), selected = c("TP53"), multiple = TRUE,
                                                       options = list(maxItems = 4))
                                 ),
                                 column(3, fileInput("custom_gb",label = HTML("OR upload your own GenBank file<b style='display:inline; color:red;' id ='q1'> [?!]</b>"),
                                                     multiple=F,accept=c('.gb','.gbk'),width = '100%')
                                 )
                             ),
                        # )
                        HTML("currently selected references"),
                        uiOutput("alignTo_new")
                    )
                )
				,column(12,tags$div(title="",HTML(paste("<br><hr><div style='display:inline;'><h4 style='display:inline;'>samples</h4> &nbsp <b style='display:inline;' id ='samples_ui'> [?]</b></div>")))
				)
				,fluidRow(
                    column(6,
                        fileInput("browser_files","upload .abi or .ab1 files with unique names",multiple=T,accept=c('.abi','.ab1'),width = '100%')),
                    column(2,
                        tags$div(title="The list of confirmed variants can be exported and saved in the form of an Excel table.",
                            downloadButton('export_btn','export variants [?]',class = "exp_btn")))
                    ,DT::dataTableOutput('samples_table')
		        )
            ),
			tabPanel('variants', value = 'main', icon = icon("search"), # http://fontawesome.io/icons/
				fluidRow(
					column(1,
					    # this should be replaced by direct interaction with graph or data table
                        tags$div(title="",
                            textInput("choose_call_pos",HTML("<div>call position <b id = 'call_pos_help'>[?]</b></div>"))),
					    tags$div(title="",
    					    selectInput("change_user_sample",HTML("<div>change variants <b id = 'change_user_help'>[?]</b></div>"),choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    selectInput("change_user_mut",NULL,choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    # textInput("change_variant","change variant to"),
                            actionButton("change_btn","change", icon = icon("exchange"), width='100%'))
					),
					column(3,
                        tags$div(title="",
    				        HTML(paste(HTML("<div>general infobox <b id = 'general_infobox_help'>[?]</b></div>"), sep=""))),
					    verbatimTextOutput("infobox")
					),
					column(2,
                        tags$div(title=""
    					    ,HTML(paste("<div>hetero indels infobox <b id = 'hetero_indel_help'>[?]</b></div>"))
    					    ,verbatimTextOutput("hetero_indel_pid")
    					)
					    ,conditionalPanel(condition = "!output.indels_present",
					        HTML("<font color=lightgrey><i>if indels detected, checkbox will appear</i></font>")
					    )
                        ,conditionalPanel(condition = "output.indels_present",
                            tags$div(title="if there are indel events above, use them to try and correct the variant calling",
                                checkboxInput("incorporate_checkbox","use detected hetero indels [?]", value = F)
                            )
					    )
					),
					column(1,
                        tags$div(title="",
                            sliderInput("mut_min",HTML("Detection limit <b id= mut_min_help> [?]</b>"), ticks=FALSE, min = 0, max = 50, value = 10, step = 0.5, round = 1)
                        )
					),
					column(1,
					       tags$div(title="",
					                sliderInput("s2n_min",HTML("min S/N <b id = min_sn_help>[?]</b>"), ticks=FALSE, min = 0, max = 10, value = 2, step = 0.1, round = 1)
					       )
					    ,sliderInput("qual_thres_to_call","min quality", ticks=FALSE, min = 0, max = 50, value = 0)
                        #sliderInput("qual_thres_to_trim","[qual thres to trim]", ticks=FALSE, min = 0, max = 60, value = 0)
					),
                    column(2
                        ,column(12,
                            tags$div(title="",HTML("ignore sequence positions <b id = trim_help >[?]</b>")))
                        ,column(6
                            ,numericInput("trim_fwd_start","fwd, up to",value = 0,min = 0,max=1000)
                            ,numericInput("trim_rev_start","rev, up to",value = 0,min = 0,max=1000)
                        )
                        ,column(6
                            ,numericInput("trim_fwd_end",  "fwd, after",value = 1000,min = 1,max=1000)
                            ,numericInput("trim_rev_end",  "rev, after",value = 1000,min = 1,max=1000)
                        )
                    ),
					column(1
					    ,checkboxInput("show_qual_checkbox","show quality", value = F)
					    ,checkboxInput("show_calls_checkbox","show calls", value = F)
					    ,conditionalPanel(condition = "output.reverse",
                                          checkboxInput("join_traces_checkbox","join traces", value = F))
                        #sliderInput("opacity_fwd","fwd trace opacity", ticks=FALSE, min = 0, max = 100, value = 100, step = 5),
                    ),
					column(1
                        ,sliderInput("max_y_p","peak height", ticks=FALSE, min = 1, max = 400, value = 100, step = 10)
					    ,conditionalPanel(condition = "output.reverse && input.join_traces_checkbox" ,
					                      sliderInput("opacity","R <trace opacity> F", ticks=FALSE, min = -100, max = 100, value = 0, step = 10))
					    ,conditionalPanel(condition = "!input.join_traces_checkbox", HTML("<div id='spacer'>  </div>"))
					    #,uiOutput('helpButton')
					    #,actionButton("toggle_help",icon = icon("question"), "hide help",class = "show_hide_help")
					    ,tags$style(HTML("#spacer{margin-top:43px;}"))
					)
				),
				fluidRow(
				    #  === following lines left out of text bc of space ===
				    # "| horizontal grey line = full sequence |"
				    #  | <font color=brown>brown dots</font> = intensity anomalies (indels?)

				    #input.toggle_help  % 2 == 0

					conditionalPanel("input.toggle_help % 2 == 0",column(12,wellPanel(tags$div(HTML(paste("<div style=\"font-family:'Inconsolata';font-size:1.1em;\">
                        <b>minimap&nbsp&nbsp&nbsp&nbsp&nbsp</b>: <font color=#4682B4>blue box = resize/move for navigation</font> | boxes = exons/introns | horizontal grey line = full sequence | verticals = variants, ref>pri>sec | <font color=red>red dotted lines</font> = filtered noisy beginnings</br>
                        <b>chromatogram</b>: click text to print info^ | crosshair + click = drag | zoom-in for extra info (coords..) | sequences, from top = ref, call/pri, mut/sec | pink verticals = variants | grey bars = codons</br>
                        <b>variants&nbsp&nbsp&nbsp&nbsp</b>: 'goto' = go to variant on chromatogram | 'x' = ignore for the session | 'confirm' = keep for the session (even if you change parameters) and make them exportable from 'samples' panel
                        </div>
  		            "), sep = "")))))
				),
				chromatographyOutput("plot"),
				DT::dataTableOutput("chosen_variants_table"),
				#conditionalPanel(condition=" output.chosen_variants_table ",downloadButton("export_btn","export")),
				br()
			),
            tabPanel('help', value = 'help', icon = icon("question"),
                         # selectInput("img_help","choose a topic:",choices=c("File Upload","Select ref from list","Load custom GenBank ref.","Scroll","Zoom","Export detected variant"))
                         selectInput("img_help","choose a topic:",choices=c("select ref from list","load custom GenBank ref","file upload","scroll","zoom","export detected variants"))
                         ,uiOutput("upload_file")
   			)
            ,
			tabPanel('hetero alignment', value = 'aln', icon = icon("sliders"),
		         verbatimTextOutput("aln"),
		         plotOutput('het_histogram',height = 600)
			)
# 			tabPanel('variant annotation', value = 'annotation', icon = icon("user-md")
# 			),
#            ,
# 			tabPanel('calls', value = 'call_table', icon = icon("table"),
# 				 shiny::dataTableOutput("call_table")
#  			)
# 			,tabPanel('intensities fwd', value = 'intens_table', icon = icon("table"),
# 		         shiny::dataTableOutput("intens_table")
# 			)
#             ,tabPanel('intensities rev', value = 'intens_table_rev', icon = icon("table"),
#                  shiny::dataTableOutput("intens_table_rev")
#             )
		),
        bsModal("modalnew", "Change name", "BUTnew", size = "small",
            HTML("Do you want to change the name?"),
            actionButton("BUTyes", "Yes"),
            actionButton("BUTno", "No")
        )
    )
)