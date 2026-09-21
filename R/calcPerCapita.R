#' Calculate per-capita variables based on GDP|PPP
#'
#' @param dt Aggregated MIF variables.
#'
#' @returns Per-capita variables.
#' @import data.table
#' @noRd

calcPerCapita <- function(dt) {

  variable <- value <- population <- conversionFactor <- newUnit <- unit <- NULL

  passengerVariables <- c(
    "ES|Transport|Pass",
    "ES|Transport|Pass with bunkers",
    "ES|Transport|Pass|Aviation",
    "ES|Transport|Bunkers|Pass|International Aviation",
    "ES|Transport|Pass|Domestic Aviation",
    "ES|Transport|Pass|Road|Bus",
    "ES|Transport|Pass|Non-motorized|Walk",
    "ES|Transport|Pass|Non-motorized|Cycle",
    "ES|Transport|Pass|Non-motorized",
    "ES|Transport|Pass|Rail",
    "ES|Transport|Pass|Rail|non-HSR",
    "ES|Transport|Pass|Rail|HSR",
    "ES|Transport|Pass|Road|LDV",
    "ES|Transport|Pass|Road",
    "ES|Transport|Pass|Road|LDV|Four Wheelers",
    "ES|Transport|Pass|Road|LDV|Three Wheelers",
    "ES|Transport|Pass|Road|LDV|Two Wheelers",
    "ES|Transport|Pass|Road|LDV|BEV",
    "ES|Transport|Pass|Road|LDV|FCEV",
    "ES|Transport|Pass|Road|LDV|Gases",
    "ES|Transport|Pass|Road|LDV|Hybrid electric",
    "ES|Transport|Pass|Road|LDV|Liquids"
  )
  freightVariables <- c(
    "ES|Transport|Freight",
    "ES|Transport|Freight with bunkers",
    "ES|Transport|Bunkers|Freight|International Shipping",
    "ES|Transport|Freight|Road",
    "ES|Transport|Freight|Domestic Shipping",
    "ES|Transport|Freight|Rail",
    "ES|Transport|Freight|Road|BEV",
    "ES|Transport|Freight|Road|FCEV",
    "ES|Transport|Freight|Road|Gases",
    "ES|Transport|Freight|Road|Liquids"
  )

  #build a lookup table with one row per variable that should get a p.c. twin
  perCapitaVariables <- data.table(
    variable = c("GDP|PPP", passengerVariables, freightVariables,
                 "Stock|Transport|Pass|Road|LDV|Four Wheelers"),
    newUnit = c("kUS$2017", rep("pkm/yr", length(passengerVariables)),
                rep("tkm/yr", length(freightVariables)),
                "cars per 1000 people"),
    conversionFactor = c(1e6, rep(1e9, length(passengerVariables) +
                                    length(freightVariables) + 1))
  )

  byColumns <- intersect(c("model", "scenario", "region", "period"), names(dt))
  populationData <- copy(dt[variable == "Population",
                            c(byColumns, "value"),
                            with = FALSE])
  populationData[, population := value * 1e6][, value := NULL]

  perCapita <- merge(dt, perCapitaVariables, by = "variable")
  perCapita <- merge(perCapita, populationData, by = byColumns)
  perCapita[, value := value / population * conversionFactor]
  perCapita[, `:=`(variable = paste0(variable, " pCap"), unit = newUnit)]
  perCapita[, c("newUnit", "conversionFactor", "population") := NULL]
  setcolorder(perCapita, names(dt))

  return(perCapita)
}
