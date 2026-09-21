test_that("configured sales technology shares are calculated", {
  testData <- new.env(parent = emptyenv())
  load(test_path("testVariableAggregation.RData"), envir = testData)

  sales <- copy(testData$vars[variable == "ES"])
  sales[, variable := "Sales"]
  aggregatedSales <- aggregateVariables(sales, testData$mapAggregation)
  shareGroups <- salesTechnologyShareGroups()
  shares <- calcShares(aggregatedSales, shareGroups)

  expectedVariables <- paste0(
    "Share|",
    unlist(Map(paste, names(shareGroups), shareGroups, MoreArgs = list(sep = "|")),
           use.names = FALSE)
  )
  expect_setequal(unique(shares$variable), expectedVariables)
  expect_length(expectedVariables, 17)
  expect_false(anyNA(shares))

  denominator <- "Sales|Transport|Pass|Road|Bus"
  busVariables <- paste(denominator, shareGroups[[denominator]], sep = "|")
  expectedShares <- merge(
    aggregatedSales[variable %in% busVariables],
    aggregatedSales[variable == denominator, .(region, period, total = value)],
    by = c("region", "period")
  )
  expectedShares[, value := value / total]
  expectedShares[, `:=`(variable = paste0("Share|", variable), unit = "-")]
  expectedShares[, total := NULL]
  actualShares <- copy(shares[variable %in% expectedShares$variable])
  setkey(actualShares, NULL)
  setkey(expectedShares, NULL)

  expect_equal(actualShares, expectedShares,
               ignore.col.order = TRUE, ignore.row.order = TRUE)
})

test_that("technology sales must add up to total sales", {
  denominator <- "Sales|Transport|Pass|Road|Bus"
  sales <- data.table(region = "World",
                      variable = c(denominator,
                                   paste(denominator, c("BEV", "Liquids"), sep = "|")),
                      unit = "million vehicles/yr",
                      period = 2030,
                      value = c(10, 4, 5))

  expect_error(calcShares(sales,
                          setNames(list(c("BEV", "Liquids")), denominator)),
               "do not sum")
})
