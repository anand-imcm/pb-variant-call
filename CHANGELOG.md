# Changelog

All notable changes to this project will be documented in this file.
## [1.2.4] - 2022-02-19

### Added
- VEP annotation for structural variants.
- VEP annotated SV VCF to TSV -> Useful in filtering/sorting annotation data.
- Github Actions upgrade

## [1.2.3] - 2022-01-31

### Added
- VEP annotated VCF to TSV -> Useful in filtering/sorting annotation data.

## [1.2.2] - 2022-01-30

### Fixed
- Coverage plot artifacts -> Sort coverage data by position

## [1.2.1] - 2022-01-26

### Fixed
- Summary module -> Improved file handling by checking if the input files exist and are non-empty before opening them.
- Summary module -> Added proper closing of file handles after reading the files.

## [1.2.0] - 2022-01-22

### Added

- New depth of coverage plot with target genomic intervals.
- Added a new input param to the workflow `region_to_plot` which contains all the genomic coordinates of MTX1, GBAP1, MTX1P1 and GBA1.

## [1.1.0] - 2023-12-19

### Added

- Added an additional summary of variants by applying a variant allele fractions (VAF) filter of greater than 0.5 to improve the confidence of the variants.

## [1.0.0] - 2023-12-16

### Added

- Initial release adapted from https://github.com/anand-imcm/wf-pb-amp