#'Report basic variable set needed to report REMIND input data and detailed transport data
#'
#' @param data List that contains at least the model results to report the basic variable set
#' @param timeResReporting Timesteps to be reported
#'
#' @returns list of intensive and extensive output variables
#' @author Johanna Hoppe
#' @import data.table
#' @export

toolReportBaseVarSet <- function(data, timeResReporting) {

  subsectorL3 <- variable <- period <- NULL

  # aggregate costs---------------------------------------------------------------------------
  aggregatedCosts <- toolReportAggregatedCosts(data$combinedCAPEXandOPEX)
  aggregatedCosts <- merge(aggregatedCosts, data$helpers$decisionTree,
                           by = intersect(names(aggregatedCosts), names(data$helpers$decisionTree)))

  # Move from sales to fleet reporting for affected variables---------------------------------
  # (in the variables named fleet other modes are still included)
  # Energy service demand on fleet level deviates from the sales level
  # regarding the share that each technology gets
  fleetESdemand <- rbind(data$ESdemandFVsalesLevel[!grepl("Bus.*|.*4W|.*freight_road.*", subsectorL3)],
                         data$fleetSizeAndComposition$fleetESdemand)
  # Energy intensity and Capital costs are tied to the construction year and have to be recalculated
  # to reflect the value for each year referring to the vehicle stock
  fleetVariables <- list(fleetEnergyIntensity = data$enIntensity,
                         fleetCapCosts = aggregatedCosts[variable == "Capital costs sales"])
  fleetData <- lapply(fleetVariables, toolReportFleetVariables,
                      data$fleetSizeAndComposition$fleetVehNumbersConstrYears, data$helpers)
  fleetCost <- rbind(fleetData$fleetCapCosts, aggregatedCosts[!variable == "Capital costs sales"])

  # Calculate final energy---------------------------------------------------------------------
  fleetFEdemand <- toolReportFE(fleetEnergyIntensity = fleetData$fleetEnergyIntensity, fleetESdemand = fleetESdemand,
                                loadFactor = data$loadFactor, hybridElecShare = data$hybridElecShare, helpers = data$helpers)

  # Split extensive and intensive variables ---------------------------------------------------
  outputVarsExt <- list(fleetESdemand = fleetESdemand,
                        fleetFEdemand = fleetFEdemand)
  outputVarsInt <- list(fleetEnergyIntensity = fleetData$fleetEnergyIntensity,
                        fleetCost = fleetCost)
  outputVars <- list(ext = outputVarsExt,
                     int = outputVarsInt)

  return(outputVars)

}
