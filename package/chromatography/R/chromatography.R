#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export

chromatography <- function(intens, intens_rev = NULL, intrexdat, calls, choices, new_sample, width = NULL, height = NULL) {

  #data = fromJSON(file=Data)
    x <- list(
        intens     = intens,
        intens_rev = intens_rev,
        intrexdat  = intrexdat,
        calls      = calls,
        choices    = choices,
        new_sample = new_sample
    )
  # create widget
    htmlwidgets::createWidget(
        name = 'chromatography',
        x,
        width = width,
        height = height,
        package = 'chromatography'
    )
    #create input binding
}
#' Widget output function for use in Shiny
#'
#' @export
chromatographyOutput <- function(outputId, width = '100%', height = '500px'){
    shinyWidgetOutput(outputId, 'chromatography', width, height, package = 'chromatography')
}
#' Widget render function for use in Shiny
#'
#' @export
renderChromatography <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) { expr <- substitute(expr) } # force quoted
    shinyRenderWidget(expr, chromatographyOutput, env, quoted = TRUE)
}
