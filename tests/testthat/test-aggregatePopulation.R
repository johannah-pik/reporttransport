test_that("population aggregation accepts unmapped passthrough regions", {
  population <- data.table(
    region = c("A", "B", "C"),
    period = 2030,
    variable = "Population",
    unit = "million",
    value = c(10, 20, 30)
  )
  regSubsetMap <- data.table(
    region = c("A", "B"),
    aggrReg = "AB"
  )

  result <- aggregatePopulation(population, regSubsetMap)

  expect_setequal(result$region, c("A", "B", "C", "AB", "World"))
  expect_equal(result[region == "AB", value], 30)
  expect_equal(result[region == "World", value], 60)
})

test_that("population aggregation skips an incompatible regional resolution", {
  population <- data.table(
    region = c("AB", "C"),
    period = 2030,
    variable = "Population",
    unit = "million",
    value = c(30, 30)
  )
  regSubsetMap <- data.table(
    region = c("A", "B"),
    aggrReg = "AB"
  )

  result <- aggregatePopulation(population, regSubsetMap)

  expect_setequal(result$region, c("AB", "C", "World"))
  expect_equal(result[region == "World", value], 60)
})
