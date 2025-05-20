## ----echo = FALSE, results = "hide"-------------------------------------------
knitr::opts_chunk$set(
  error = FALSE,
  fig.width = 7,
  fig.height = 5)
set.seed(1)

## ----echo = FALSE, results = "asis"-------------------------------------------
writeLines(c(
  "```c",
  readLines(system.file("include/cinterpolate/cinterpolate.h",
                        package = "cinterpolate")),
  "```"))

## ----echo = FALSE, results = "asis"-------------------------------------------
writeLines(c(
  "```plain",
  readLines(system.file("example/DESCRIPTION",
                        package = "cinterpolate")),
  "```"))

## ----echo = FALSE, results = "asis"-------------------------------------------
writeLines(c(
  "```plain",
  readLines(system.file("example/NAMESPACE",
                        package = "cinterpolate")),
  "```"))

## ----echo = FALSE, results = "asis"-------------------------------------------
writeLines(c(
  "```c",
  readLines(system.file("example/src/testing.c",
                        package = "cinterpolate")),
  "```"))

