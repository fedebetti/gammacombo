# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.1.0/)
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.4.0 [Unreleased]

### Added

### Changed
* Removed stateless classes `ColorBuilder`, `FitResultDump` and `TGraphTools`
* Moved `float` -> `double` for all variables except for the ones stored in
  `TTree`s and the ones that are `float` in ROOT (`Minuit` internally uses
  `double`).

### Fixed
* Const correctness of (most) class methods.
* Delete copy constructors and copy assignment operators for resource-managing
  classes.
* Headers are now self-sufficient.
