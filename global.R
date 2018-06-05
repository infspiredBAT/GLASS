library(data.table)

general_g_calls                 <- NULL             #annotated basecall data
general_g_intens                <- NULL             #intensities file
general_g_intens_rev            <- NULL             #optional reverse intensities file
general_g_single_rev            <- NULL
general_g_intrexdat             <- NULL             #intrex data used in graphs
general_g_glassed_ref           <- NULL
general_g_glassed_cod           <- NULL
general_g_custom_ref            <- NULL
general_g_custom_cod            <- NULL
general_g_glassed_snp           <- NULL
general_g_alignTo_options       <- c("TP53","ATM","NOTCH1","CALR")
general_g_alignTo_description   <- list("TP53"   = "<b>TP53</b> description"
                                        , "ATM"    = "<b>ATM</b> description"
                                        , "NOTCH1" = "<b>NOTCH1</b> description"
                                        , "CALR"   = "<b>CALR</b> description")
general_g_choices               <- NULL
general_g_noisy_neighbors       <- NULL
general_g_view                  <- NULL
general_g_selected              <- NULL
general_g_selected_goto_index   <- 0

general_g_calls                 <- NULL             #annotated basecall data
general_g_intens                <- NULL             #intensities file
general_g_intens_rev            <- NULL             #optional reverse intensities file
general_g_single_rev            <- NULL
general_g_intrexdat             <- NULL             #intrex data used in graphs
general_g_glassed_ref           <- NULL
general_g_glassed_cod           <- NULL
general_g_custom_ref            <- NULL
general_g_custom_cod            <- NULL
general_g_glassed_snp           <- NULL
general_g_alignTo_options       <- c("TP53","ATM","NOTCH1","CALR")
general_g_alignTo_description   <- list("TP53"   = "<b>TP53</b> description"
                                        , "ATM"    = "<b>ATM</b> description"
                                        , "NOTCH1" = "<b>NOTCH1</b> description"
                                        , "CALR"   = "<b>CALR</b> description")
general_g_choices               <- NULL
general_g_noisy_neighbors       <- NULL
general_g_view                  <- NULL
general_g_selected              <- NULL
general_g_selected_goto_index   <- 0
general_g_max_y                 <- NULL
general_g_hetero_calls          <- 0
general_g_hetero_indel_pid      <- 0
general_g_hetero_ins_tab        <- NULL
general_g_hetero_del_tab        <- NULL
general_g_expected_het_indel    <- NULL
general_g_minor_het_insertions  <- NULL
general_load_id                 <- NULL
general_g_stored_het_indels     <- list()
general_files_info              <- ""
general_g_indels_present        <- FALSE
general_g_qual_present          <- FALSE
general_g_brush_fwd             <- 0
general_g_brush_rev             <- 0
general_g_not_loaded            <- ""
general_g_reactval              <- reactiveValues()
general_g_reactval$updateVar    <- 0
general_g_refs_avail            <- c("-","TP53","NOTCH1","ATM","CALR","Custom")
general_g_files         <- data.table(FWD_name=character(), FWD_file=character(), REV_name=character(), REV_file=character(),
                                       REF=character(), mut_min=numeric(), qual_thres_to_call=numeric(), s2n_min=numeric(),
                                       show_calls_checkbox=logical(), join_traces_checkbox=logical(), max_y_p=numeric(), 
                                       opacity=numeric(), incorporate_checkbox=logical(), loaded=logical(), calls=character(),
                                       status=character(), id=integer(), brush_fwd_start=numeric(), brush_fwd_end=numeric(), 
                                       brush_rev_start=numeric(), brush_rev_end=numeric(), coding=character(), protein=character(), 
                                       VAF=character(), dbSNP=character(), dbSNP_id=character()) #g_files = represents the data table containing individual samples their info

if (("dev" %in% list.files())) {

    source("procAbi.R", local = TRUE)
    source("helpers.R", local = TRUE)
    source("samples.R", local = TRUE)

    general_global_vars <- ls()
    for (x in general_global_vars[grepl("general_", general_global_vars)]) {
      assign(gsub("general_", "", x), get(x))
    }
}
