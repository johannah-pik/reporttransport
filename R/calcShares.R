#' Calculate shares based on aggregated variables ("...4W|BEV" with "...|4W" as the denominator)
#'
#' @param dt  containing reporting variables.
#' @param shareGroups Named list of denominator variables and their categories.
#' @param tolerance Numerical tolerance for the share check.
#'
#' @returns dt containing the calculated share variables.
#' @import data.table
#' @noRd

calcShares <- function(dt, shareGroups, tolerance = 1e-6) {

  variable <- value <- total <- numeratorTotal <- denominatorTotal <- NULL
  idColumns <- setdiff(names(dt), c("variable", "unit", "value"))

  shares <- Map(function(denominatorVariable, technologies) {
    numeratorVariables <- paste(denominatorVariable, technologies, sep = "|")
    expectedVariables <- c(denominatorVariable, numeratorVariables)
    missingVariables <- setdiff(expectedVariables, unique(dt$variable))
    if (length(missingVariables) > 0) {
      stop(paste("Missing variables:", paste(missingVariables, collapse = ", ")))
    }

    totals <- copy(dt[variable == denominatorVariable,
                      c(idColumns, "value"), with = FALSE])
    setnames(totals, "value", "total")
    shareRows <- merge(dt[variable %chin% numeratorVariables], totals,
                       by = idColumns)

    check <- shareRows[,
                       list(numeratorTotal = sum(value),
                            denominatorTotal = unique(total)),
                       by = idColumns]
    difference <- abs(check$numeratorTotal - check$denominatorTotal)
    if (any(difference > tolerance * pmax(1, abs(check$denominatorTotal)))) {
      stop(paste("Technology sales do not sum to", denominatorVariable))
    }

    shareRows[, value := ifelse(total == 0, 0, value / total)]
    shareRows[, `:=`(variable = paste0("Share|", variable), unit = "-")]
    return(shareRows[, names(dt), with = FALSE])
  }, names(shareGroups), shareGroups)

  return(rbindlist(shares))
}

salesTechnologyShareGroups <- function() {
  return(list(
    "Sales|Transport|Pass|Road|LDV|Two Wheelers" = c("BEV", "Liquids"),
    "Sales|Transport|Pass|Road|LDV|Three Wheelers" = c("BEV", "Liquids"),
    "Sales|Transport|Pass|Road|Bus" = c("BEV", "FCEV", "Gases", "Liquids"),
    "Sales|Transport|Pass|Road|LDV|Four Wheelers" =
      c("BEV", "FCEV", "Gases", "Hybrid electric", "Liquids"),
    "Sales|Transport|Freight|Road" = c("BEV", "FCEV", "Gases", "Liquids")
  ))
}
