# GLASS
### assisted and standardised assessment of gene variations from Sanger data

Live version available at http://bat.infspire.org/genomepd/glass

Requirements
--------------

* RStudio or <href link="https://www.rstudio.com/products/shiny/shiny-server/">Shiny Server</href>

Tested on :
<verbatim>
``` R
sessionInfo()
R version 3.3.3 (2017-03-06)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X El Capitan 10.11.6

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] seqinr_3.3-6              zoo_1.8-0                 xlsx_0.5.7                xlsxjars_0.6.1           
 [5] rJava_0.9-8               sangerseqR_1.10.0         shinyBS_0.61              DT_0.2                   
 [9] stringi_1.1.5             stringr_1.2.0             chromatography_0.0.0.9000 htmlwidgets_0.8          
[13] rjson_0.2.15              data.table_1.10.4         shinyjs_0.9               shiny_1.0.3              
[17] Biostrings_2.42.1         XVector_0.14.1            IRanges_2.8.2             S4Vectors_0.12.2         
[21] BiocGenerics_0.20.0      

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.10    magrittr_1.5    zlibbioc_1.20.0 lattice_0.20-35 xtable_1.8-2    R6_2.2.1        tools_3.3.3    
 [8] grid_3.3.3      miniUI_0.1.1    htmltools_0.3.6 ade4_1.7-6      yaml_2.1.14     digest_0.6.12   mime_0.5       
[15] jsonlite_1.4    httpuv_1.3.3
```
