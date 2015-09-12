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

		HTML("<font size=3em><b>genomePD/glass</b></font><font size=2.5em> | <i>dev</i> | Pal Bystry Reigl Krejci Demko and Darzentas @ <a href=http://bat.infspire.org>bat.infspire.org</a></font><font size=1.0em> | CESNET/MetaCentrum | CEITEC MU</font>"),

		hr(),

		tabsetPanel(id = 'tab',
			tabPanel('variants', value = 'main', icon = icon("search"), # http://fontawesome.io/icons/
				fluidRow(
					column(2,
					    fileInput("select_file","",multiple=T,accept=c('.abi','.ab1'))
					),
					column(1,
					    # this should be replaced by direct interaction with graph or data table
					    textInput("choose_variance","type seq pos"),
					    selectInput("change_peak","user_mod it to",choices=c("A","T","C","G"),selected="",selectize=F,size=1),
					    actionButton("execute_btn","change", icon = icon("exchange"))
					),
					column(3,
					    HTML(paste("messages and info will appear here", sep="")),
					    verbatimTextOutput("variance_info")
					),
					column(2,
					    # sliderInput("rm7qual_thres","set rm7qual thres for trimming", ticks=FALSE, min = 0, max = 50, value = 12),
					    sliderInput("qual_thres","set qual thres for low qual", ticks=FALSE, min = 0, max = 50, value = 10)
					    # sliderInput("aln_min","set min coverage for alignment", ticks=FALSE, min = 0, max = 1, value = 0.2)
					),

					column(2,
					    sliderInput("max_y_p","set intens peak height", ticks=FALSE, min = 0, max = 200, value = 100),
              sliderInput("intens_guide_line","set intens guideline rel height", ticks=FALSE, min = 0, max = 100, value = 50)
					),
          column(2,
              sliderInput("opacity_fwd","set forward trace opacity", ticks=FALSE, min = 0, max = 100, value = 100),
              sliderInput("opacity_rev","set reverse trace opacity", ticks=FALSE, min = 0, max = 100, value = 100)
          )
				),

				chromatographyOutput("plot"),

				shiny::dataTableOutput("chosen_variances_table"),

				conditionalPanel(condition=" output.chosen_variances_table ",
	                 downloadButton("export_btn","export"),
	                 actionButton("delete_btn","delete", icon = icon("times"))
				)
	 		),
			tabPanel('annotation', value = 'annotation', icon = icon("user-md")
			),
			tabPanel('calls', value = 'call_table', icon = icon("table"),
				 shiny::dataTableOutput("call_table")
 			),
			tabPanel('intensities fwd', value = 'intens_table', icon = icon("table"),
		         shiny::dataTableOutput("intens_table")
			),
            tabPanel('intensities rev', value = 'intens_table_r', icon = icon("table"),
                 shiny::dataTableOutput("intens_table_r")
            )
		)

    ))
