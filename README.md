# MicrobiomeExperiment

<!-- badges: start -->

[![R-CMD-check-Bioc-devel](https://github.com/FelixErnst/MicrobiomeExperiment/workflows/R-CMD-check-Bioc-devel/badge.svg)](https://github.com/FelixErnst/MicrobiomeExperiment/actions)
[![Codecov test
coverage](https://codecov.io/gh/FelixErnst/MicrobiomeExperiment/branch/master/graph/badge.svg)](https://codecov.io/gh/FelixErnst/MicrobiomeExperiment?branch=master)

<!-- badges: end -->

This project is aimed to provide a `SummarizedExperiment` infrastructure
for microbiome experiment data. It is based on the `TreeSummarizedExperiment`
package and adds additional functionality from the microbiome analysis world to
the context of working with `SummarizedExperiment`.

The `MicrobiomeExperiment` package currently serves as a wrapper for the project
pulling in:

- classes from the [`TreeSummarizedExperiment`](https://github.com/fionarhuang/TreeSummarizedExperiment) package
- methods for data wrangling from the [`mia`](https://github.com/FelixErnst/mia)  package
- methods for microbiome data visualization from the [`miaViz`](https://github.com/microbiome/miaViz) package

Workflows, which span the whole packages, should end up in `MicrobiomeExperiment`
as individual vignettes.

# Contribution

Feel free to contribute to the individual packages accordingly or in case of 
example workflows to `MicrobiomeExperiment`

## Technical aspects

Let's use a git flow kind of approach. Development version should be done 
against the `master` branch and then merged to `release` for release. 
(https://guides.github.com/introduction/flow/)

# Code of conduct

Please note that the MicrobiomeExperiment project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
