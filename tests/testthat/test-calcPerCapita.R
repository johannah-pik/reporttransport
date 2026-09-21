test_that("per-capita variables are calculated from population", {
  dt <- data.table(
    region = "World",
    variable = c("Population", "GDP|PPP", "ES|Transport|Pass|Road|Bus",
                 "Stock|Transport|Pass|Road|LDV|Four Wheelers"),
    unit = c("million", "billion US$2017", "billion pkm/yr", "million vehicles"),
    period = 2030,
    value = c(2, 100, 4, 1)
  )

  perCapita <- calcPerCapita(dt)
  expected <- data.table(
    region = "World",
    variable = c("ES|Transport|Pass|Road|Bus pCap", "GDP|PPP pCap",
                 "Stock|Transport|Pass|Road|LDV|Four Wheelers pCap"),
    unit = c("pkm/yr", "kUS$2017", "cars per 1000 people"),
    period = 2030,
    value = c(2000, 50, 500)
  )
  setkey(perCapita, NULL)

  expect_equal(perCapita, expected,
               ignore.col.order = TRUE, ignore.row.order = TRUE)
})
