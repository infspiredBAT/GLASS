#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export
#'


chromatography <- function(new_sample,
                           show_calls = FALSE,
                           show_qual = FALSE,
                           filters = TRUE,
                           num_samples,
                           samples,
                           width = NULL, height = '600px') {

    print("new sample")
    print(new_sample)

    x <- list(
        new_sample  = new_sample,
        resize      = FALSE,
        show_calls  = show_calls,
        show_qual   = show_qual,
        filters     = filters,
        num_samples = num_samples,
        samples     = samples
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
chromatographyOutput <- function(outputId, width = '100%', height = '700px'){
    shinyWidgetOutput(outputId, 'chromatography', width, height, package = 'chromatography')
}
#' Widget render function for use in Shiny
#'
#' @export
renderChromatography <- function(expr, env = parent.frame(), quoted = FALSE) {
    if (!quoted) { expr <- substitute(expr) } # force quoted
    shinyRenderWidget(expr, chromatographyOutput, env, quoted = TRUE)
}
