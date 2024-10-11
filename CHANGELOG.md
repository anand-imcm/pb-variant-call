# Changelog

All notable changes to this project will be documented in this file.

## [1.2.7] - 2024-10-11

### Added

- Script to generate a summary table for depth of coverage per gene.

### Updated

- Docker registry from DockerHub to GitHub to avoid DockerHub's rate limit restrictions.
- Dynamic runtime and bootdisk for handling `PAPI error code 9`.

## [1.2.6] - 2024-06-24

### Fixed

- Input to accept both fastq.gz and .fastq file format to the workflow.

## [1.2.5] - 2024-02-19

### Added

- Variant Filter based on QUAL > 30
- Filtered variants summary module now generates a variant-level summary .tsv file for all the variants based on the new filter criteria.
- VAF summary removed.

## [1.2.4] - 2024-02-19

### Added

- VEP annotation for structural variants.
- VEP annotated SV VCF to TSV -> Useful in filtering/sorting annotation data.
- Github Actions upgrade

## [1.2.3] - 2024-01-31

### Added

- VEP annotated VCF to TSV -> Useful in filtering/sorting annotation data.

## [1.2.2] - 2024-01-30

### Fixed

- Coverage plot artifacts -> Sort coverage data by position

## [1.2.1] - 2024-01-26

### Fixed

- Summary module -> Improved file handling by checking if the input files exist and are non-empty before opening them.
- Summary module -> Added proper closing of file handles after reading the files.

## [1.2.0] - 2024-01-22

### Added

- New depth of coverage plot with target genomic intervals.
- Added a new input param to the workflow `region_to_plot` which contains all the genomic coordinates of MTX1, GBAP1, MTX1P1 and GBA1.

## [1.1.0] - 2023-12-19

### Added

- Added an additional summary of variants by applying a variant allele fractions (VAF) filter of greater than 0.5 to improve the confidence of the variants.

## [1.0.0] - 2023-12-16

### Added

- Initial release adapted from https://github.com/anand-imcm/wf-pb-amp