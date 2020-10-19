# MicrobiomeExperiment

This project is aimed to provide a `SummarizedExperiment` infrastructure
for microbiome experiment data. It is based on the `TreeSummarizedExperiment`
package and adds additional functionality from the microbiome analysis world to
the context of working with `SummarizedExperiment`.

Currently this includes:

- recognition of taxonomic information in `rowData`
- row and column merging
- agglomeration by taxonomic information (like phyloseq)

# ToDo

In this package

- ~addition of reference sequences slots~
- addition of specialized taxonomic container, if necessary
- ...

In other packages

- alpha and beta diversity measure (Mia)
- plotting functions
- ...

# Contribution

Feel free to contribute.

## Technical aspects

Let's use a git flow kind of approach. Development version should be done 
against the `master` branch and then merged to `release` for release. 
(https://guides.github.com/introduction/flow/)

# Code of conduct

Please note that the MicrobiomeExperiment project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
