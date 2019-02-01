# timeseriesNEM
A method for mapping a non-interventional time series onto a static nested effects model, as inferred using the Bioconductor package `nem`.

## Installing

Please make sure that you have the package `devtools` installed. Then run:
`devtools::install_github("cbg-ethz/timeseriesNEM")`

### Creating the vignette

In order to install the vignette, add the option `build_vignettes = TRUE` to the above command. If the vignette is still not created, the following command should do the trick:
`remotes::install_github("cbg-ethz/timeseriesNEM", build_opts = c("--no-resave-data", "--no-manual"), force = TRUE)`
