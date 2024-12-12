# SpQ - QSpike Tools reinvented - electrophysiology extracellular multichannel batch and parallel preprocessor
#    Copyright (C) 2024 Michele GIUGLIANO <michele.giugliano@unimore.it> and contributors.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>

function parse_toml_files(OUTPUT)::settings

   #-- TOML --------------------------------------
   config = joinpath(OUTPUT, "config.toml")

   if !isfile(config)
       @error "$config not found."
       exit(1)
   end

   configfile  = TOML.parsefile(config)
   LFP    = configfile["lfp"]["save_lfp"] == "y"
   fmax   = configfile["lfp"]["fmax"]
   rate   = configfile["lfp"]["rate"]

   detect = configfile["detection"]["detect"] == "y"
   fmin_d = configfile["detection"]["fmin"]
   fmax_d = configfile["detection"]["fmax"]
   stdmin = configfile["detection"]["stdmin"]
   stdmax = configfile["detection"]["stdmax"]
   event  = configfile["detection"]["type"] == "both" ? 0 : configfile["detection"]["type"] == "pos" ? 1 : -1
   ref    = configfile["detection"]["refractory"]          # ms

   shapes = configfile["sorting"]["save_shapes"] == "n"
   fmin_s = configfile["sorting"]["fmin"]
   fmax_s = configfile["sorting"]["fmax"]
   dpre   = configfile["sorting"]["delta_pre"]             # ms
   dpost  = configfile["sorting"]["delta_post"]            # ms
   factor = configfile["sorting"]["interpolation_factor"]
   spline = configfile["sorting"]["cubic_spline"] == "y"
   #----------------------------------------------

   meta = joinpath(OUTPUT, "meta.toml")
   if !isfile(meta)
       @error "$meta not found. Something must have gone wrong with the bash script."
       exit(1)
   end

   infofile = TOML.parsefile(meta)
   Nchans   = infofile["info"]["Nchans"]     # Number of channels
   Nsamples = infofile["info"]["Nsamples"]   # Number of samples per channel

   Tick     = infofile["chans"]["Tick"]      # Sampling interval in us
   ADZero   = infofile["chans"]["ADZero"]    # ADZero
   ConversionFactor = infofile["chans"]["ConversionFactor"] # D/A conversion factor
   Exponent = infofile["chans"]["Exponent"]  # Exponent for the conversion factor
   #----------------------------------------------
   srate = 1E6 / Tick                          # Sampling rate in Hz
   c     = ConversionFactor * 10. ^ Exponent;  # Conversion from AD to physical units (V)
   d     = - ADZero * c;                       # Conversion from AD to physical units (V)
   #----------------------------------------------

   # Convert fmin_d, fmax_d, fmin_s, fmax_s, dpre, dpost to Float32, before calling the functions
   bpfilt   = SpQ.prepare_bandpass(Float32(fmin_d), Float32(fmax_d), Float32(srate));
   bpfilt_s = SpQ.prepare_bandpass(Float32(fmin_s), Float32(fmax_s), Float32(srate));
   lpfilt   = SpQ.prepare_lowpass(Float32(fmax), Float32(srate));

   #-- Settings struct ---------------------------
  return SpQ.settings(OUTPUT, LFP, fmax, rate, detect, fmin_d, fmax_d, stdmin, stdmax, event, ref, shapes, fmin_s, fmax_s, dpre, dpost, factor, spline, Nchans, Nsamples, Tick, ADZero, ConversionFactor, Exponent, srate, c, d, lpfilt, bpfilt, bpfilt_s);

end # function parse_toml_files ----------------


