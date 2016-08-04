library(shiny)
#library(shinyjs)
library(data.table)
library(rjson)
library(htmlwidgets)
library(chromatography)
library(stringr)
library(stringi)
library(DT)
#source("JS.R")

shinyUI(
    fluidPage(
#        useShinyjs(),
#        extendShinyjs(text = jsDeleteRow),
#        extendShinyjs(text = jsSwapRow),
        theme = "simplex3.3.6.css", # http://bootswatch.com/ | sandstone/simplex/flatly/darkly
		tags$head(
		    includeCSS("www/samples.css"),
		    tags$head(HTML("<link href='https://fonts.googleapis.com/css?family=Inconsolata' rel='stylesheet' type='text/css'>")),
            tags$title("GLASS"),
		    #tags$head(tags$script(src="selectize.min.js")),
		    #tags$script(incljs),
            # script for input file cleaning
		    tags$script('
		                window.onbeforeunload = function() {
		                return "You will need to start over if you navigate away!";
		                };
		                '),
            tags$script('
              		Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {
              		var el = $("#" + x);
              		el.replaceWith(el = el.clone(true));
              		var id = "#" + x + "_progress";
              		$(id).css("visibility", "hidden");
              		});
              		'),
            tags$style(HTML(".DTFC_LeftBodyLiner { width: 100% !important; }"))
		),

		       #conditionalPanel(condition = "output.show_sample_brows==TRUE",
		       #          absolutePanel(id = "samples", class = "panel panel-default", fixed = TRUE,
		        #      draggable = TRUE, top = 200, left = 130, right = "auto", bottom = "auto",
		        #      width = 1100, height = "auto",
		        #
		        #      h2("Sample Browser"),
		        #
		        #      DT::dataTableOutput('samples_table'),
		        #      fileInput("browser_files",NULL,multiple=T,accept=c('.abi','.ab1'))
		        #         )
		        #),
		HTML(paste("&nbsp&nbsp<b><font size=6em>GLASS</font></b><font size=3em> | <a href=http://bat.infspire.org target=_blank>bat.infspire.org</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.ceitec.eu/ceitec-mu/medical-genomics/rg34 target=_blank>Medical Genomics Group @ CEITEC MU</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.ericll.org target=_blank>European Research Initiative on CLL / ERIC</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.igcll.org target=_blank>IgCLL group</a></font> | CESNET/MetaCentrum")),
		fluidRow(
			#column(1, selectInput("gene_of_interest",NULL,choices=list("ATM"="ATM","NOTCH1"="NOTCH1","TP53"="TP53"),selected="TP53",multiple=FALSE,selectize=F,size=1)),
			#column(1, selectInput("gene_of_interest",NULL,choices=list("TP53"="TP53"),selected="TP53",multiple=FALSE,selectize=F,size=1)),
            #column(1,
            #    conditionalPanel(condition = "input.gene_of_interest == 'TP53'",
            #        actionButton("ex_btn","example",icon = icon("play"),class="btn btn-info",style="width:100%;height:20px;padding:0;margin-top:8px;"))),
			#column(1,actionButton("ex_btn","example",icon = icon("play"),class="btn btn-info",style="width:100%;height:20px;padding:0;margin-top:8px;margin-bottom:4px;")),
            #column(3,wellPanel(fluidRow(
            #    column(5,
    		#	       tags$div(title="please make sure either or both (in case of paired i.e. forward and reverse) files have an \"F\" or \"R\" before the .abi/.ab1 file extension,\ne.g. my_sampleF.abi and/or my_sampleR.abi\n\nuse \"*R.abi\" even if loading a single reverse file!",
					# column(12,wellPanel(tags$div(HTML(paste("<div style=\"font-family:'Inconsolata';font-size:1em;\">
    		#	                tags$div(HTML(paste("<div style=\"font-family:'Inconsolata';font-size:1em;\">up to two ABI files,</br>name fwd and rev as:</br>*<strong style=\"color: red;\">F</strong>.abi *<strong style=\"color: red;\">R</strong>.abi [?]</div>"), sep = "")))),
    		#	column(7, fileInput("select_file",NULL,multiple=T,accept=c('.abi','.ab1')))
			#))),
			#column(2,actionButton("mng_samples_btn","Manage/Load samples",style="width:100%;height:20px;padding:0;margin-top:8px;")),
			column(2,  HTML(paste("&nbsp<i>assisted and standardised assessment</br>&nbsp&nbsp&nbspof gene aberrations from Sanger data</i></br></br>"))),
			column(9, htmlOutput("files")),
			column(1,HTML(paste(' <a href="javascript:void(0)" id="app-disclaimer-link" onclick="$(\'#disclaimer-modal\').modal(\'show\')">version 0.2.2 <br>(2016-Aug-4)</a>
                                <!-- Update log -->
                                <div class="modal fade" id="disclaimer-modal" tabindex="-1" role="dialog">
			                    <div id="disclaimer-modal-content" class="modal-dialog" role="document">
			                    <div class="modal-content">
			                    <div class="modal-header" style="padding-top:10px;">
			                    <h4>update log</h4>
			                    </div>
			                    <div class="modal-body" style="padding-top:0px;padding-bottom:0px;font-size:12px;color:rgb(132,132,132)">
                                <b>ver 0.2.2 (2016-Aug-4)</b>
                                <ul>
                                <li>Fixed bug causing single strand samples to crash. (Result of previous update.)</li>
                                <li>Small change in navigation tool opacity.</li>
                                </ul>
                                <b>ver 0.2.1 (2016-Aug-2)</b>
			                    <ul>
			                    <li>Brush (navigation tool) changed to steal blue color to represent glass (magnifier) and added shadow to make it "pop-out".</li>
			                    <li>Variant indicators are now strand specific.</li>
                                <li>Variant indicator scales to peak width.</li>
			                    </ul>
                                <b>ver 0.2.0 (2016-Jul-29)</b>
			                    <ul>
			                    <li>Mainly brush (navigation tool) appearence and behaviour.</li>
			                    <li>More seamless transition.</li>
			                    <li>Box represenging zoomed in area.</li>
                                <li>Filters of beginnings of reads stylised in chromatogram and also minimap.</li>
			                    </ul>
                                <b>ver 0.1.1 (2016-Jul-27)</b>
			                    <ul>
			                    <li>Filters of beginnings and ends of reads visualised in minimap.</li>
                                <li>Reference "case sensitive"; exons shown in uppercase, everything else is in lowercase.</li>
                                <li>It\'s now possible to load samples without a reference.</li>
			                    </ul>
			                    <b>ver 0.1.0 (2016-Jul-06)</b>
			                    <ul>
                                <li>Implemented brush-able filtering-out of variants from the beginning of reads (the filter is at the end of the reverse sequences since these are reverse complemented). </li>
                                <li>Added reference selector, with selected references used to autodetect strandedness - not using all references speeds up the upload process.</li>
                                <li>Added dialog preventing accidental "navigate away".</li>
			                    <li>Added CALR reference.</li>
			                    <li>Started this log.</li>
			                    </ul>
	   		                    </div>
			                    <div class="modal-footer" style="clear-both">
			                    <button type="button" class="btn btn-default icon-button-sm btn btn-default btn-raised" data-dismiss="modal">Close</button>
			                    </div>
			                    </div>
			                    </div>
			                    </div>
			                    ')))
		),
		tabsetPanel(id = 'tabs',
            tabPanel('samples',value = 'smpl_brws',icon = icon("flask"),
				fluidRow(
					column(12,wellPanel(tags$div(HTML(paste("<div style=\"font-family:'Inconsolata';font-size:1.1em;\">
                        <b>(very) quick guide</b>: (0) hovering over '[?]'s will provide quick help tips</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(1) upload ABI files from supported references, forward and/or reverse, and with unique names</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(2) after processing and auto-detection of most file properties, pair / unpair / swap / delete / change reference as necessary</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(3) press blue 'load' button</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(4) when done in 'variants' panel (find separate instructions there), confirmed variants will appear under status and are exportable with green 'export' button
                        </div>"), sep = ""))))),
                fluidRow(
                    column(12,tags$div(title="Please make sure the files have unique names. Uploading a file with the same name as one of the files in the table will be ignored.",
                        HTML(paste("<b>upload and process ABI files</b> [?]"))))
                ),
                fluidRow(
                    column(12,tags$div(title="Uploaded ABI files are aligned against the references of the selected genes and their orientation is automatically detected. Any of the available references can be selected. The more references are selected the longer the upload process wil take.",
                        checkboxGroupInput("alignTo", "select reference(s) to autodetect [?]", c("TP53","ATM","NOTCH1","CALR"), selected = c("TP53","CALR"),inline=TRUE)))
                ),
                fluidRow(
                    column(6,
                        fileInput("browser_files",NULL,multiple=T,accept=c('.abi','.ab1'),width = '100%')),
                    column(2,
                        tags$div(title="The list of confirmed variants can be exported and saved in the form of an Excel table.",
                            downloadButton('export_btn','export variants [?]',class = "exp_btn")))
		        ),
		        DT::dataTableOutput('samples_table')
            ),
			tabPanel('variants', value = 'main', icon = icon("search"), # http://fontawesome.io/icons/
				fluidRow(
					column(1,
					    # this should be replaced by direct interaction with graph or data table
                        tags$div(title="the absolute position of a call, nothing to do with genomic or codon numbering\n\nyou can either type a number, or it will show by interacting with glass",
                            textInput("choose_call_pos","call position [?]")),
					    # selectInput("change_peak","user_mod it to",choices=list(empty="","-","A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),selected="",selectize=F,size=1),
                        tags$div(title="change the 1st/2nd or major/minor or sample/mutation variants to...",
    					    selectInput("change_user_sample","change variants [?]",choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    selectInput("change_user_mut",NULL,choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    # textInput("change_variant","change variant to"),
                            actionButton("change_btn","change", icon = icon("exchange"), width='100%'))
					),
					column(4,
                        tags$div(title="information about the current position, also highlighted in light blue, including all nucleotides (called, reference, mutated), coordinates, qualities (Q), peak %s, signal-to-noise (S/N) ratios, etc",
    				        HTML(paste("general infobox [?]", sep=""))),
					    verbatimTextOutput("infobox")
					),
					column(2,
                        tags$div(title="whether we expect heterozygous indels in the data,\nthe % identity of the alignment between the primary and secondary sequences (get this as high as you can),\nand insertion / deletion (in that order) counts\n\nus expecting an indel does not mean there is one: try any suggested peak%, keep %id high, and use your best judgement after studying the chromatogram\n\nif indels are detected, a checkbox to use them will appear under this infobox",
    					    HTML(paste("hetero indels infobox [?]")),verbatimTextOutput("hetero_indel_pid")
    					),
#                         tags$div(title="distinct insertion / deletion (in that order) events, and their lengths in nt",
#     					    HTML(paste("hetero ins/dels [?]")),verbatimTextOutput("hetero_indel_tab")
#     					),
					    conditionalPanel(condition = "!output.indels_present",
					        HTML("<font color=lightgrey><i>if indels detected, checkbox will appear</i></font>")
					    ),
					    conditionalPanel(condition = "output.indels_present",
                            tags$div(title="if there are indel events above, use them to try and correct the variant calling",
                                checkboxInput("incorporate_checkbox","use detected hetero indels [?]", value = F)
                            )
					    )
					),
					column(1,
                        tags$div(title="min % of peak for mutation to be called",
                            sliderInput("mut_min","mut: min peak% [?]", ticks=FALSE, min = 0, max = 50, value = 20, step = 0.5, round = 1)
                        ),
                        tags$div(title="min signal to noise ratio for mutation to be called",
    					    sliderInput("s2n_min","mut: min S/N [?]", ticks=FALSE, min = 0, max = 10, value = 2, step = 0.1, round = 1)
                        )
					),
					column(1
					    ,sliderInput("qual_thres_to_call","min quality", ticks=FALSE, min = 0, max = 50, value = 0)
                        #sliderInput("qual_thres_to_trim","[qual thres to trim]", ticks=FALSE, min = 0, max = 60, value = 0)
					),
                    column(1),
					column(1
					    ,checkboxInput("show_calls_checkbox","   show calls", value = F)
					    ,conditionalPanel(condition = "output.reverse",
                                          checkboxInput("join_traces_checkbox","join traces", value = F))
                        #sliderInput("opacity_fwd","fwd trace opacity", ticks=FALSE, min = 0, max = 100, value = 100, step = 5),
                    ),
					column(1
                        ,sliderInput("max_y_p","peak height", ticks=FALSE, min = 1, max = 200, value = 100, step = 10)
					    ,conditionalPanel(condition = "output.reverse && input.join_traces_checkbox" ,
					                      sliderInput("opacity","R <trace opacity> F", ticks=FALSE, min = -100, max = 100, value = 0, step = 10))
					)
				),
				fluidRow(
				    # "| horizontal grey line = full sequence |" left out of text bc of space
					column(12,wellPanel(tags$div(HTML(paste("<div style=\"font-family:'Inconsolata';font-size:1.1em;\">
                        <b>minimap&nbsp&nbsp&nbsp&nbsp&nbsp</b>: <font color=#4682B4>blue box = resize/move for navigation</font> | boxes = exons/introns | verticals = variants, ref>pri>sec | <font color=red>red dotted lines</font> = filter out noisy beginnings | <font color=brown>brown dots</font> = intensity anomalies (indels?)</br> 
                        <b>chromatogram</b>: click text to highlight call and show info^ | sequences from top = reference, call/primary, mutation/secondary | striped verticals = indicators e.g. variants | grey bars = quality</br>
                        <b>variants&nbsp&nbsp&nbsp&nbsp</b>: 'goto' = go to variant on chromatogram | 'x' = ignore for the session | 'confirm' = keep for the session (even if you change parameters) and make them exportable from 'samples' panel
                        </div>
  		            "), sep = ""))))),
				chromatographyOutput("plot"),
				DT::dataTableOutput("chosen_variants_table"),
				#conditionalPanel(condition=" output.chosen_variants_table ",downloadButton("export_btn","export")),
				br()
			)
#            ,
#			tabPanel('hetero alignment', value = 'aln', icon = icon("sliders"),
#		         verbatimTextOutput("aln"),
#		         plotOutput('het_histogram',height = 600)
#			)
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
		)
    )
)
