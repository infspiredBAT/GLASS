#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
chromatography <- function(ins, helperdat, width = NULL, height = NULL) {
  #data = fromJSON(file=Data)
    x <- list(
        data = ins,
        helperdat = helperdat
    )
  # create widget
    htmlwidgets::createWidget(
        name = 'chromatography',
        x,
        width = width,
        height = height,
        package = 'chromatography'
    )
}
#' Widget output function for use in Shiny
#'
#' @export
chromatographyOutput <- function(outputId, width = '100%', height = '400px'){
    shinyWidgetOutput(outputId, 'chromatography', width, height, package = 'chromatography')
}
#' Widget render function for use in Shiny
#'
#' @export
renderChromatography <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) { expr <- substitute(expr) } # force quoted
    shinyRenderWidget(expr, chromatographyOutput, env, quoted = TRUE)
}