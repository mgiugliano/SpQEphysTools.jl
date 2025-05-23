# SpQEphysTools - QSpike Tools reinvented - electrophysiology extracellular multichannel batch and parallel preprocessor
#    Copyright (C) 2024 Michele GIUGLIANO <michele.giugliano@unimore.it> and contributors.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>

const ZPG = DSP.Filters.ZeroPoleGain{:z,ComplexF32,ComplexF32,Float32}; # Type alias

@with_kw mutable struct settings
  OUTPUT::String = ""  # Output directory
  LFP::Bool = false  # Extract LFP
  type::String = "" # Type of data (MCS or 3BRAIN)
  fmax::Float32 = 0.0 # Low-pass cutoff frequency for LFP
  rate::Float32 = 0.0 # Decimation rate for LFP

  detect::Bool = false # Detect spikes
  fmin_d::Float32 = 0.0 # Low-pass cutoff frequency for spike detection
  fmax_d::Float32 = 0.0 # High-pass cutoff frequency for spike detection
  stdmin::Float32 = 0.0 # Threshold for spike detection
  stdmax::Float32 = 0.0 # Threshold for spike detection
  event::Int = 0  # Event type for spike detection (1: positive, -1: negative, 0: both)
  ref::Float32 = 0.0 # Refractory period for spike detection (ms)

  shapes::Bool = false  # Extract spike shapes
  fmin_s::Float32 = 0.0  # Low-pass cutoff frequency for spike shapes
  fmax_s::Float32 = 0.0  # High-pass cutoff frequency for spike shapes
  dpre::Float32 = 0.0  # Pre window for spike shapes (ms)
  dpost::Float32 = 0.0  # Post window for spike shapes (ms)
  factor::Int = 0   # Factor for spike shapes
  spline::Bool = false # Spline interpolation for spike shapes

  Nchans::Int = 0 # Number of channels
  Nsamples::Int = 0 # Number of samples

  Tick::Float32 = 0.0 # Time of the first sample
  ADZero::Float32 = 0.0 # Zero level of the A/D converter
  ConversionFactor::Float32 = 0.0 # Conversion factor (D/A)
  Exponent::Float32 = 0.0         # Exponent for the conversion factor
  srate::Float32 = 0.0         # Sampling rate
  c::Float32 = 0.0         # Conversion factor (D/A)
  d::Float32 = 0.0         # Conversion factor (D/A)
  #n_work::Int
  lpfilt::ZPG = ZPG([0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], 1.0)
  bpfilt::ZPG = ZPG([0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], 1.0)
  bpfilt_s::ZPG = ZPG([0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], 1.0)
end  # struct settings ----------------------------


