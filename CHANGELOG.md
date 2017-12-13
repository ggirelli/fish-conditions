# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [unreleased]
### Changed
- Renamed to MUD-FISH.

### Fixed
- First hybridization duplex type now passed to 2nd hybridization instead of being fixed to DNA:RNA.

## [2.0.0]
### Fixed
- Updated oligo_melting library to version with fixed thermodynamic table labeling issue.
- Temperature score function now takes into account different interaction type for main target and additional ones.
- Temperature score function now allows different concentration of additional targets.
- Temperature score function now takes into account different probe sizes (in oligo).

### Changed
- Moved settings confirmation code to common library.
- Auxiliary output moved to `aux` folder.
- Single picked table, for convenience.
- Single temperature condition score table, fore convenience.
- Increased the number of decimals for normalized score.

### Added
- Plot score distributions.
- Saving command line also for multi-probe script.
- Confirm output directory overwrite before running.
- Now logging the best normalized score.

## [1.2.0]
### Fixed
- MacOS support for sed versions.
- MacOS and Windows support for different file formats.

### Added
- Now multiple probes with same color sequence are supported.

## [1.1.1]
### Fixed
- Parallelization support on Mac.

## [1.1.0]
### Added
- Support for Mac in the documentation.
- Support for Magnesium concentration effect.
- Normalized score.
- Updated oligo-melting with automatic R package installation.
- Added more info to picked.tsv table.

## [1.0.0]

* [unreleased] https://github.com/ggirelli/mud-fish
* [2.0.0] https://github.com/ggirelli/mud-fish/releases/tag/v2.0.0
* [1.2.0] https://github.com/ggirelli/mud-fish/releases/tag/v1.2.0
* [1.1.1] https://github.com/ggirelli/mud-fish/releases/tag/v1.1.1
* [1.1.0] https://github.com/ggirelli/mud-fish/releases/tag/v1.1.0
* [1.0.0] https://github.com/ggirelli/mud-fish/releases/tag/v1.0.0
