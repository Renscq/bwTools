test_that("bwTools declares intentional data.table semantics", {
  expect_true(isTRUE(get(".datatable.aware", envir = asNamespace("bwTools"), inherits = FALSE)))
})
