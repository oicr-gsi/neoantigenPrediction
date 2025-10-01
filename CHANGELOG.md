# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-09-30
### Changed
- Changed wdl to use callable scripts from module neopipe/1.1.0 from inline code.

## [1.0.2] - 2025-08-11
### Added
- New input parameter vep_buffer_size (default: 500) to control the VEP buffer size, passed via --vep_buffer_size.

### Fixed
- HLA parsing logic now excludes invalid alleles (".") by adding checks for allele1 != "." and allele2 != "."

## [1.0.1] - 2025-06-04
### Changed
- HLACalls now stores a single file and caller instead of arrays.

## [1.0.0] - 2025-03-24
### Added
- [GRD-843](https://jira.oicr.on.ca/browse/GRD-843), Initial verion of the wdl along with README and vidarr files
