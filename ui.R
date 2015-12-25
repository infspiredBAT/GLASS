library(shiny)
library(data.table)
library(rjson)
library(htmlwidgets)
library(chromatography)
library(stringr)
library(stringi)
#require(DT)

shinyUI(

    fluidPage(theme = "simplex.css", # http://bootswatch.com/ | sandstone/simplex/flatly/darkly
		tags$head(
            tags$title("genomePD/glass"),
            # script for input file cleaning
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

		HTML(paste("&nbsp&nbsp<b><font size=4em>genomePD/ </font><font size=6em>glass</font></b><font size=3em> | <a href=http://bat.infspire.org target=_blank>bat.infspire.org</a> &nbsp.&nbsp <a href=http://www.ceitec.eu/ceitec-mu/medical-genomics/rg34 target=_blank>the Medical Genomics Group @ CEITEC MU</a> &nbsp.&nbsp <a href=http://www.ericll.org target=_blank>the European Research Initiative on CLL / ERIC</a> &nbsp.&nbsp <a href=http://www.igcll.org target=_blank>the IgCLL group</a></font> | CESNET/MetaCentrum")),
		fluidRow(
    		br(),
			column(1, selectInput("gene_of_interest","",choices=list("ATM"="ATM","NOTCH1"="NOTCH1","TP53"="TP53"),selected="TP53",multiple=FALSE,selectize=F,size=1)),
            column(1,
                conditionalPanel(condition = "input.gene_of_interest == 'TP53'",
                    actionButton("ex_btn","example",icon = icon("play"),class="btn btn-info",style="width:100%;height:20px;padding:0;margin-top:8px;"))),
			column(2, fileInput("select_file","",multiple=T,accept=c('.abi','.ab1'))),
			column(8, verbatimTextOutput("files"))
		),
		tabsetPanel(id = 'tab',
			tabPanel('variants', value = 'main', icon = icon("search"), # http://fontawesome.io/icons/
				fluidRow(
					column(1,
					    # this should be replaced by direct interaction with graph or data table
                        tags$div(title="the absolute position of a call, nothing to do with genomic or codon numbering\n\nyou can either type a number, or it will show by interacting with glass",
                            textInput("choose_call_pos","call position [?]")),
					    # selectInput("change_peak","user_mod it to",choices=list(empty="","-","A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),selected="",selectize=F,size=1),
                        tags$div(title="change the 1st/2nd or major/minor or sample/mutation variants to...",
    					    selectInput("change_user_sample","change variants [?]",choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    selectInput("change_user_mut","",choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
    					    # textInput("change_variant","change variant to"),
                            actionButton("change_btn","change", icon = icon("exchange"), width='100%'))
					),
					column(4,
				        HTML(paste("general infobox", sep="")),
					    verbatimTextOutput("infobox")
					),
					column(2,
                        tags$div(title="whether we can detect heterozygous indels in the data,\nthe % identity of the alignment between the sequences of the primary and mutation calls,\nand distinct insertion / deletion (in that order) events, and their lengths in nt\n\nif this low and there are many variants reported, there might be a heterozygous indel",
    					    HTML(paste("hetero indels infobox [?]")),verbatimTextOutput("hetero_indel_pid")
    					),
#                         tags$div(title="distinct insertion / deletion (in that order) events, and their lengths in nt",
#     					    HTML(paste("hetero ins/dels [?]")),verbatimTextOutput("hetero_indel_tab")
#     					),
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
					    ,sliderInput("qual_thres_to_call","min quality", ticks=FALSE, min = 0, max = 50, value = 20)
                        #sliderInput("qual_thres_to_trim","[qual thres to trim]", ticks=FALSE, min = 0, max = 60, value = 0)
					),
                    column(1),
					column(1
					    ,checkboxInput("show_calls_checkbox","show calls", value = F)
					    ,conditionalPanel(condition = "output.reverse",
                                          checkboxInput("join_traces_checkbox","join traces", value = F))
                        #sliderInput("opacity_fwd","fwd trace opacity", ticks=FALSE, min = 0, max = 100, value = 100, step = 5),
                    ),
					column(1
                        ,sliderInput("max_y_p","peak height", ticks=FALSE, min = 1, max = 200, value = 100, step = 10)
					    ,conditionalPanel(condition = "output.reverse",
					                      sliderInput("opacity","R <trace opacity> F", ticks=FALSE, min = -100, max = 100, value = 0, step = 10))
					)
				),
				br(),

				chromatographyOutput("plot"),
				shiny::dataTableOutput("chosen_variants_table"),
				conditionalPanel(condition=" output.chosen_variants_table ",
				    downloadButton("export_btn","export")
				),
				br()
			),
			tabPanel('hetero alignment', value = 'aln', icon = icon("sliders"),
		         verbatimTextOutput("aln"),
		         plotOutput('het_histogram',height = 600)
			),
# 			tabPanel('variant annotation', value = 'annotation', icon = icon("user-md")
# 			),
			tabPanel('calls', value = 'call_table', icon = icon("table"),
				 shiny::dataTableOutput("call_table")
 			)
# 			,tabPanel('intensities fwd', value = 'intens_table', icon = icon("table"),
# 		         shiny::dataTableOutput("intens_table")
# 			)
#             ,tabPanel('intensities rev', value = 'intens_table_rev', icon = icon("table"),
#                  shiny::dataTableOutput("intens_table_rev")
#             )
		)
    )
)
