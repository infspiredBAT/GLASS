#' <Add Title>
#'
#' <Add Description>
#'
#' @import htmlwidgets
#'
#' @export

chromatography <- function(intens, intens_rev = NULL,single_rev, intrexdat, calls, choices, new_sample, noisy_neighbors, show_calls = FALSE,show_qual = FALSE,qual_present,brush_fwd_start,brush_fwd_end,brush_rev_start,brush_rev_end, width = NULL, height = NULL) {

  #data = fromJSON(file=Data)
    x <- list(
        intens     = intens,
        intens_rev = intens_rev,
        single_rev = single_rev,
        intrexdat  = intrexdat,
        calls      = calls,
        choices    = choices,
        new_sample = new_sample,
        resize     = FALSE,
        noisy_neighbors = noisy_neighbors,
        show_calls      = show_calls,
        show_qual       = show_qual,
        qual_present    = qual_present,
        brush_fwd_start = brush_fwd_start,
        brush_fwd_end   = brush_fwd_end,
        brush_rev_start = brush_rev_start,
        brush_rev_end   = brush_rev_end
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
