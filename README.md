# SpQEphysTools.jl - Main library for SpQ Project

This library is developed as the backbone of the data analysis tool set [SpiQ](https://github.com/mgiugliano/SpiQ).
It contains a set of routines to perform initial preprocessing and analysis of electrophysiological recordings,
acquired extracellularly by substrate-integrated microelectrode arrays (MEAs).

It implements state-of-the-art algorithms and techniques in peak-detection (see the works of Prof. Dr. Juan Quiroga) and
extend a previous approach (i.e. [QSpikeTools](https://doi.org/10.3389/fninf.2014.00026)) originally based on MATLAB and
aimed at batch parallel processing of data on HPC.

This new Julia version is polished and minimal, and builds on top of Julia's excellent performance and simplicity.

Please refer to [SpiQ](https://github.com/mgiugliano/SpiQ) project for more information. Collaborators and pull requests are
welcome!

Prof. Dr. Ing. Michele GIUGLIANO ([www](https://www.giugliano.info)).



[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mgiugliano.github.io/SpQEphysTools.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mgiugliano.github.io/SpQEphysTools.jl/dev/)
[![Build Status](https://github.com/mgiugliano/SpQEphysTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mgiugliano/SpQEphysTools.jl/actions/workflows/CI.yml?query=branch%3Amain)
