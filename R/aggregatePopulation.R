#' Aggregate population for reporting
#'
#' Population should be able to handle different regional resolutions from (H12 and EU21).
#' Aggregate detailed regions covered by the supplied map,
#' and calculate World from the unaggregated population to avoid double counting.
#'
#' @param population Pop data in long format.
#' @param regSubsetMap Mapping with `region` and `aggrReg` columns.
#'
#' @returns Pop data containing the original, aggregated, and World regions.
#' @noRd

aggregatePopulation <- function(population, regSubsetMap) {
  region <- value <- NULL

  populationBase <- copy(population)
  populationWorld <- populationBase[, .(value = sum(value)),
                                    by = setdiff(names(populationBase), c("region", "value"))]
  populationWorld[, region := "World"]

  mapRegions <- unique(regSubsetMap$region)
  populationSubset <- populationBase[region %chin% mapRegions]

  populationAggregated <- NULL
  if (length(mapRegions) > 0L &&
        setequal(unique(populationSubset$region), mapRegions)) {
    populationAggregated <- as.data.table(
      aggregate_map(populationSubset, regSubsetMap, by = "region")
    )
  }

  return(rbindlist(
    list(populationBase, populationAggregated, populationWorld),
    use.names = TRUE,
    fill = TRUE
  ))
}
