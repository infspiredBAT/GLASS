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

		HTML("<font size=3em><b>genomePD/glass</b></font><font size=2.5em> | <i>dev</i> | Pal Bystry Reigl Krejci Demko Malcikova and Darzentas @ <a href=http://bat.infspire.org>bat.infspire.org</a></font><font size=1.0em> | CESNET/MetaCentrum | CEITEC MU</font>"),

		hr(),

		tabsetPanel(id = 'tab',
			tabPanel('variants', value = 'main', icon = icon("search"), # http://fontawesome.io/icons/
				fluidRow(
					column(2,
					    fileInput("select_file","",multiple=T,accept=c('.abi','.ab1')),
					    verbatimTextOutput("files")
					),
					column(1,
					    # this should be replaced by direct interaction with graph or data table
					    textInput("choose_variance","call position"),
					    #selectInput("change_peak","user_mod it to",choices=list(empty="","-","A","T","G","C","S","W","R","Y","K","M","B","V","H","D","N"),selected="",selectize=F,size=1),
					    selectInput("change_user_sample","sample variant to",choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
					    selectInput("change_user_mut","mutant variant to",choices=list("",deletion="-","A","T","G","C","S (G or C)"="S","W (A or T)"="W","R (A or G)"="R","Y (C or T)"="Y","K (G or T)"="K","M (A or C)"="M","B (C or G or T)"="B","V (A or C or G)"="V","H (A or C or T)"="H","D (A or G or T)"="D","N"),selected="",selectize=F,size=1),
                        actionButton("change_btn","change", icon = icon("exchange"))
					),
					column(4,
				        HTML(paste("messages and info will appear here", sep="")),
					    verbatimTextOutput("infobox")
					),
					column(1,
                        # HTML(paste("hetero calls:")),verbatimTextOutput("hetero_calls"),
					    HTML(paste("hetero aln %id:")),verbatimTextOutput("hetero_indel_pid"),
					    HTML(paste("hetero ins/dels:")),verbatimTextOutput("hetero_indel_tab"),
                        checkboxInput("incorporate_checkbox","incorporate", value = T)
					),
					column(1,
                        sliderInput("mut_min","min peak% for mut", ticks=FALSE, min = 0, max = 50, value = 20, step = 0.5, round = 1),
					    sliderInput("s2n_min","min signal/noise", ticks=FALSE, min = 0, max = 10, value = 2, step = 0.1, round = 1)
					),
					column(1,
                        sliderInput("qual_thres_to_call","qual thres to call", ticks=FALSE, min = 0, max = 50, value = 14),
                        sliderInput("qual_thres_to_trim","[qual thres to trim]", ticks=FALSE, min = 0, max = 50, value = 12)
					),
					column(1
                        #sliderInput("opacity_fwd","fwd trace opacity", ticks=FALSE, min = 0, max = 100, value = 100, step = 5),
                    ),
					column(1,
                        sliderInput("opacity","trace opacity", ticks=FALSE, min = -100, max = 100, value = 0, step = 10),
                        sliderInput("max_y_p","rel peak height", ticks=FALSE, min = 0, max = 200, value = 100, step = 10)
					)
				),
				br(),

				chromatographyOutput("plot"),

				shiny::dataTableOutput("chosen_variances_table"),

				conditionalPanel(condition=" output.chosen_variances_table ",
				    actionButton("reset_btn","reset to ref", icon = icon("times"))
				    ,downloadButton("export_btn","export")
                    ,actionButton("confirm","confirm prediction",icon = icon("check"))
				),

				br()
			),
			tabPanel('hetero alignment', value = 'aln', icon = icon("sliders"),
		         verbatimTextOutput("aln")
			),
			tabPanel('variant annotation', value = 'annotation', icon = icon("user-md")
			),
			tabPanel('calls', value = 'call_table', icon = icon("table"),
				 shiny::dataTableOutput("call_table")
 			),
			tabPanel('intensities fwd', value = 'intens_table', icon = icon("table"),
		         shiny::dataTableOutput("intens_table")
			),
            tabPanel('intensities rev', value = 'intens_table_rev', icon = icon("table"),
                 shiny::dataTableOutput("intens_table_rev")
            )
		)
    )
)
