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
            tags$style(HTML(".DTFC_LeftBodyLiner { width: 100% !important; }")),
		    HTML("<style type='text/css'>
                .form-control{
                    background-color: white !important;
                    //margin-left: 4px !important;
                }
                .btn-file{
                   
                }
                .exp_btn{
                    display:none;
                }
                .irs-bar, .irs-bar-edge{
                    visibility: hidden;
                }
                .irs-slider {
                    background-color: #777;
                }
                .zoom{
                    cursor: move;
                    fill: none;
		            pointer-events: all;
                }
                
		         </style>")
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
		fluidRow(
			column(1, HTML("&nbsp&nbsp<b><font size=6em>GLASS</font></b>")),
			column(1, HTML(paste('<div  style="padding:.3em .5em"><a href="javascript:void(0)" id="app-disclaimer-link" onclick="$(\'#disclaimer-modal\').modal(\'show\')" > 
                                ver 0.2.15 <br> (2017-Jun-7) </a></div>
                                <!-- Update log -->
                                
                                <div class="modal fade" id="disclaimer-modal" tabindex="-1" role="dialog">
			                    <div id="disclaimer-modal-content" class="modal-dialog" role="document">
			                    <div class="modal-content">
			                    <div class="modal-header" style="padding-top:10px;">
			                    <h4>update log</h4>
			                    </div>
			                    <div class="modal-body" style="padding-top:0px;padding-bottom:0px;font-size:12px;color:rgb(132,132,132)">
                                <b>ver 0.2.15 (2017-Jun-7)</b>
                                <ul>
			                     <li>New dialogues to help with custom GenBank files and stability improvement for GenBank file loading.</li>
			                     </ul>
                                <b>ver 0.2.14 (2017-Apr-11)</b>
                                <ul>
			                     <li>Handling error from empty variants table.</li>
                                  <li>Corrected trimming filter for double stranded variants.</li>
			                     </ul>
                                <b>ver 0.2.13 (2017-Apr-9)</b>
                                <ul>
			                     <li>Fixed broken status reporting in samples table.</li>
			                     </ul>
                                <b>ver 0.2.12 (2017-Apr-3)</b>
                                <ul>
			                     <li>Help images </li>
			                     <li>Custom GenBank file experimental feature introduced.</li>
                                 <li>References shown in tiles with additional information provided.</li>
                                 <li>Widgets for setting filters added.</li>
                                 <li>"show qualities" issue fixed.</li>
			                     </ul>
                                <b>ver 0.2.11 (2017-Mar-22)</b>
                                <ul>
			                     <li>Improvments in error handling. </li>
                                <li>Compatibility issue with new version of the data tables library addressed. </li>
			                    </ul>
                                <b>ver 0.2.10 (2017-Jan-13)</b>
                                <ul>
			                     <li>Redefined the formula for estimating the resolution. New unit BasePerPixel should be consistent across different input sequence lengths and window sizes. </li>
			                     </ul>
                                <b>ver 0.2.9 (2016-Nov-17)</b>
                                <ul>
			                     <li>Scrolling events ignored on zoom (zoom on scroll was not seamless).</li>
			                     <li>Reordered svg elements. Order of elements ~ z coordinate. </li>
			                     </ul>
                                <b>ver 0.2.8 (2016-Nov-3)</b>
                                <ul>
			                     <li>UI modifications. Better interactivity (Click and drag in the graph area).</li>
			                     <li>Filtered beginnings of reads now set to a fixed value. </li>
			                     </ul>
                                <b>ver 0.2.7 (2016-Oct-7)</b>
                                <ul>
			                     <li>A bug preventing the loading of some single nucleotide sequences fixed.</li>
			                     <li>Correction in the way VAF (variant allele frequency) is estimated.</li>
			                     </ul>
                                <b>ver 0.2.6 (2016-Sep-29)</b>
                                <ul>
			                    <li>Small changes in the table of samples formating for a more comprehensive representation of variants and a more unified style in both the samples table and the final exported excel table.</li>
                                <li>Added proper naming of protein variants for p.(=) and p.? cases.</li>
                                </ul>
                                <b>ver 0.2.5 (2016-Sep-10)</b>
                                <ul>
			                    <li>dbSNP annotation added for TP53 now uses exact matching instead of position matching.</li>
                                <li>Added button to show/hide the help infobox in the variants panel.</li>
                                <li>Titles, email us, table formatting, quick guide</li>
                                <li>Modified export button behaviour; hide when nothing to export.</li>
                                <li>Delete button in samples table behaviour fixed.</li>
                                <li>Enabled the delete button on the example file.</li>
			                    </ul>
                                <b>ver 0.2.4 (2016-Aug-20)</b>
                                <ul>
			                    <li>Position based dbSNP annotation added for TP53. (Column "dbSNP" in the variants table.)</li>
			                    </ul>
                                <b>ver 0.2.3 (2016-Aug-11)</b>
                                <ul>
                                <li>"No reference" upload fixed and documented.</li>
                                <li>Forward beginnings filter initial position updated.</li>
                                </ul>
                                <b>ver 0.2.2 (2016-Aug-4)</b>
                                <ul>
                                <li>Fixed bug causing single strand samples to crash (result of previous update).</li>
                                <li>Small change in minimap navigation box opacity.</li>
                                </ul>
                                <b>ver 0.2.1 (2016-Aug-2)</b>
			                    <ul>
			                    <li>Minimap navigation box coloured steel blue (~ magnifier glass) and shadow added to make it "pop-out".</li>
			                    <li>Variant indicators are now strand specific.</li>
                                <li>Variant indicator scales to peak width.</li>
			                    </ul>
                                <b>ver 0.2.0 (2016-Jul-29)</b>
			                    <ul>
			                    <li>Minimap navigation box appearence and behaviour.</li>
			                    <li>More seamless chromatogram transitions.</li>
			                    <li>Box representing zoomed-in area.</li>
                                <li>Filters of beginnings of reads stylised in chromatogram and also minimap.</li>
			                    </ul>
                                <b>ver 0.1.1 (2016-Jul-27)</b>
			                    <ul>
			                    <li>Filters of beginnings and ends of reads visualised in minimap.</li>
                                <li>Reference now case sensitive: exons in uppercase, everything else in lowercase.</li>
                                <li>It\'s now possible to load samples without a reference, just untick all genes.</li>
			                    </ul>
			                    <b>ver 0.1.0 (2016-Jul-06)</b>
			                    <ul>
                                <li>Implemented brush-able filtering-out of variants from the beginning of reads (the filter is at the end of the reverse sequences since these are reverse complemented). </li>
                                <li>Added reference selector, with selected references used to autodetect strandedness - not using all references speeds up the upload process.</li>
                                <li>Added dialog preventing accidental "navigate away".</li>
			                    <li>Added CALR reference.</li>
			                    <li>Update log started.</li>
			                    </ul>
	   		                    </div>
			                    <div class="modal-footer" style="clear-both">
			                    <button type="button" class="btn btn-default icon-button-sm btn btn-default btn-raised" data-dismiss="modal">close</button>
			                    </div>
			                    </div>
			                    </div>
			                    </div>
			                    '))),
			column(10, HTML("<br><font size=2em>&nbsp&nbsp | R&D by <a href=http://bat.infspire.org target=_blank>bat.infspire.org</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.ceitec.eu/ceitec-mu/medical-genomics/rg34 target=_blank>Medical Genomics @ CEITEC MU</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.ericll.org target=_blank>European Research Initiative on CLL / ERIC</a> &nbsp<font size=0.9em>&</font>&nbsp <a href=http://www.igcll.org target=_blank>IgCLL group</a></font>&nbsp&nbsp | IT by <a href=https://metavo.metacentrum.cz/en/ target='_blank'>CESNET/MetaCentrum</a>&nbsp&nbsp | @ <a href=mailto:bat@infspire.org target='_blank'>email us</a>"))
		),
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
			column(2, HTML(paste("&nbsp<i>assisted and standardised assessment</br>&nbsp&nbsp&nbspof gene variations from Sanger data</i></br></br>"))),
			column(9, htmlOutput("files"))
		),
		tabsetPanel(id = 'tabs',
            tabPanel('samples',value = 'smpl_brws',icon = icon("flask"),
				fluidRow(
					column(12,wellPanel(tags$div(HTML(paste("<div style=\"font-family:'Inconsolata';font-size:1.1em;\">
                        <b>quick guide</b>: (0) hovering over '[?]'s will provide quick help tips</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(1) upload ABI files from supported references, forward and/or reverse, and with unique names - if you tick relevant references beforehand, GLASS will align and auto-orientate against them</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(2) after processing and auto-detection of most file properties, pair / unpair / swap / delete / change reference as necessary</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(3) click the blue 'play' button to load the file(s)</br>&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp(4) when done in 'variants' panel (find separate instructions there), confirmed variants will appear under 'status' and are exportable with the green 'export variants' button
                        </div>"), sep = ""))))),
                fluidRow(
                    column(12,tags$div(title="",HTML(paste("<div style='display:inline;'><h3 style='display:inline;'>References </h3><b style='display:inline;' id ='ref_help'>[?]</b></div>")))
                           ,bsTooltip(id = "ref_help","Set of working references <hr> Displayed as cards are the references against which new samples are aligned (auto detect reference) and also appeare in the \"reference dropdown menu\" in the Samples table below.", placement = "right", trigger = "hover",
                                      options = NULL)
                           )
                    
                ),
				
                fluidRow(
                    column(12,tags$div(title="Uploaded ABI files are aligned against the references of the selected genes and their orientation is automatically detected. Any of the available references can be selected.

Selecting NONE of the references will assign no reference to your chromatogram. It is still possible to view it or change the reference afterwards.

The more references are selected the longer the upload process wil take."
                        #,checkboxGroupInput("alignTo", " select reference(s) to autodetect [?]", c("TP53","ATM","NOTCH1","CALR","Custom"),selected = c("TP53"),inline=TRUE)
                    ),
                        #######################
                        # alignTo_new (start) #
                        #######################
                        
                        uiOutput("alignTo_new"),
                        uiOutput("alignTo_new_debug"),
                        #####################
                        # alignTo_new (end) #
                        #####################
                        #wellPanel(
                            checkboxInput("showWP3", HTML("Manage References <i class='fa fa-wrench' aria-hidden='true'></i>")),
                            conditionalPanel(condition="input.showWP3",
                                             fluidRow(
                                                 HTML("<br>"),
                                                 column(3, 
                                                        selectizeInput("additionalRefs", " Add reference from currated list,",c("TP53","ATM","NOTCH1","CALR"), selected = c("TP53","NOTCH1"), multiple = TRUE,
                                                                       options = list(maxItems = 4))
                                                 ),
                                                 column(3, fileInput("custom_gb",label = HTML("OR upload your own GenBank file. <b style='display:inline; color:red;' id ='q1'> [?!]</b> "
                                                                                            ),
                                                                     multiple=F,accept=c('.gb','.gbk'),width = '100%'),
                                                                     bsTooltip(id = "q1","Currently an experimental feature. Help imporove it by reporting errors.<hr> Genomic coordinates are not extracted. <hr> intorn/exon numbers may be different if transcript starts with an untranslated exon.", placement = "top", trigger = "hover",
                                                                                                                                 options = NULL)
                                                       )
                                             )
                            )
                        #)
                    )
                        
                )
				,column(12,tags$div(title="",HTML(paste("<br><br><div style='display:inline;'><h3 style='display:inline;'>Samples </h3><b style='display:inline;' id ='samples_ui'>[?]</b></div><br><br>")))
				        ,bsTooltip(id = "samples_ui","Please make sure the files have unique names. Uploading a file with the same name as one of the files in the table will be ignored.", placement = "right", trigger = "hover",
				                   options = NULL)
				)
				,fluidRow(
                    column(6,
                        fileInput("browser_files",NULL,multiple=T,accept=c('.abi','.ab1'),width = '100%')),
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
                        tags$div(title="the absolute position of a call, nothing to do with genomic or codon numbering\n\nyou can either type a number, or it will show by interacting with glass",
                            textInput("choose_call_pos","call position [?]")),
					    # selectInput("change_peak","user_mod it to",choices=list(empty="","-","A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),selected="",selectize=F,size=1),
                        tags$div(title="change the 1st/2nd or major/minor or sample/mutation variants to...",
    					    selectInput("change_user_sample","change variants [?]",choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    selectInput("change_user_mut",NULL,choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    # textInput("change_variant","change variant to"),
                            actionButton("change_btn","change", icon = icon("exchange"), width='100%'))
					),
					column(3,
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
                    column(2
                        ,column(12,h4("Read trimming"))
                        ,column(6
                            ,numericInput("trim_fwd_start","start of fwd",value = 20,min = 0,max=1000)
                            ,numericInput("trim_rev_start","start of rev",value = 20,min = 0,max=1000)
                        )
                        ,column(6
                            ,numericInput("trim_fwd_end","end of fwd",value = 20,min = 0,max=1000)
                            ,numericInput("trim_rev_end","end of rev",value = 20,min = 0,max=1000)
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
                        ,sliderInput("max_y_p","peak height", ticks=FALSE, min = 1, max = 200, value = 100, step = 10)
					    ,conditionalPanel(condition = "output.reverse && input.join_traces_checkbox" ,
					                      sliderInput("opacity","R <trace opacity> F", ticks=FALSE, min = -100, max = 100, value = 0, step = 10))
					    ,conditionalPanel(condition = "!input.join_traces_checkbox", HTML("<div id='spacer'>  </div>"))
					    #,uiOutput('helpButton')
					    #,actionButton("toggle_help",icon = icon("question"), "hide help",class = "show_hide_help")
					 #   ,HTML('<button class="btn btn-default action-button show_hide_help shiny-bound-input" type="button" onclick="$(\'#disclaimer-modal2\').modal(\'show\')"><i class="fa fa-question" ></i> show help</button>
					 #           
                     #          <div class="modal fade" id="disclaimer-modal2" tabindex="-1" role="dialog">
			         #           <div id="disclaimer-modal-content" class="modal-dialog" role="document">
					 #          <div class="modal-content">
#
					 #         <div class="modal-header" style="padding-top:10px;">
					 #         <h4>Help</h4>
					 #         </div>
#
					 #         <div class="modal-body" style="padding-top:0px;padding-bottom:0px;font-size:12px;color:rgb(132,132,132)">
                     #           <br>
                     #         <ul class="nav nav-tabs">
                     #           <li class="active"><a data-toggle="tab" href="#home">Navigation</a></li>
                     #           <li><a data-toggle="tab" href="#menu1">Peaks</a></li>
                     #           <li><a data-toggle="tab" href="#menu2">Menu 2</a></li>
                     #         </ul>
#
                     #         <div class="tab-content">
                     #           <div id="home" class="tab-pane fade in active">
                     #             <h3>Scrolling along the sequence.</h3>
                     #               <img src="navbar_scroll.gif" alt="Navbar zoom" >
                     #               <h3>Zooming.</h3>
                     #               <img src="navbar_zoom.gif" alt="Navbar zoom" >
                     #           </div>
                     #           <div id="menu1" class="tab-pane fade">
                     #             <h3>Menu 1</h3>
                     #             <p>Some content in menu 1.</p>
                     #           </div>
                     #           <div id="menu2" class="tab-pane fade">
                     #             <h3>Menu 2</h3>
                     #             <p>Some content in menu 2.</p>
                     #           </div>
                     #         </div>
                     #          
                     #           </div>
					 #         <div class="modal-footer" style="clear-both">
					 #         <button type="button" class="btn btn-default icon-button-sm btn btn-default btn-raised" data-dismiss="modal">close</button>
					 #         </div>
#
					 #         </div>
					 #         </div>
					 #         </div>
					 #         ')
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
            tabPanel('Help', value = 'help', icon = icon("question"),
                         selectInput("img_help","Chose a topic:",choices=c("File Upload","Select ref from list","Load custom GenBank ref.","Scroll","Zoom","Export detected variant"))
                         ,uiOutput("upload_file")
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
		),
        bsModal("modalnew", "Change name", "BUTnew", size = "small",
        HTML("Do you want to change the name?"),
        actionButton("BUTyes", "Yes"),
        actionButton("BUTno", "No")
)
    )
)
