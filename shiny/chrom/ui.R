require(shiny)
require(data.table)
require(rjson)
require(htmlwidgets)
require(chromatography)
require(stringr)
require(stringi)
#require(DT)

shinyUI(
    fluidPage(theme = "simplex.css",

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
		HTML("<font size=3em><b>genomePD/glass</b></font><font size=2.5em> | <i>dev</i> | Pal Bystry Demko and Darzentas @ <a href=http://bat.infspire.org>bat.infspire.org</a></font><font size=1.0em> | CESNET/MetaCentrum | CEITEC MU</font>"),
		hr(),

        fluidRow(
            column(1,
                fileInput("select_file","",multiple=F,accept=c('.abi','.ab1'))
            ),
            column(2,
                # this should be replaced by direct interaction with graph or data table
                textInput("choose_variance","type sequence position to see info"),
                selectInput("change_peak","and change the call to",choices=c("A","T","C","G"),selected="",selectize=F,size=1),
                actionButton("execute_btn","change")
            ),
            column(4,
            	   # HTML(paste("info", sep="<br/>")),
            	   verbatimTextOutput("variance_info")
            ),
            column(3,
            	   sliderInput("max_y_p","set peak height",min = 0, max = 200, value = 100)
            )
        ),

        chromatographyOutput("plot"),

        # urcite bude conditional
        #        conditionalPanel(
        #            condition = " output.chosenCheckboxes == true ",
        #            DT::dataTableOutput("table2"),
        shiny::dataTableOutput("chosen_variances_table"),
        #        ),

        # vypis vybranych htmlOutput()
        downloadButton("export_btn","export"),
        #        actionButton("zoom_btn","zoom"),
        #        actionButton("edit_btn","edit"),
        actionButton("delete_btn","delete")
    )
)