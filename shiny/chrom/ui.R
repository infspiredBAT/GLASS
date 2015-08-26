require(shiny)
require(data.table)
require(rjson)
require(htmlwidgets)
require(chromatography)
require(stringr)
#require(DT)
    
shinyUI(
    fluidPage(theme = "simplex.css",
        fluidRow(
            column(6,
                # neco
                fileInput("select_file","select file",multiple=F,accept=c('.abi','.ab1'))
            ),
            column(6,
                # this should be replaced by direct interaction with graph or data table
                textInput("choose_variance","choose variance by number..."),
                verbatimTextOutput("variance_info"),
                selectInput("change_peak","and change it to",choices=c("A","T","C","G"),selected="A",selectize=F,size=1),
                actionButton("execute_btn","execute")
            )
        ),
        
        # chromatograf 
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