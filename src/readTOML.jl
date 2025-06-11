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

function parse_toml_files(OUTPUT)::settings

    #-- TOML --------------------------------------
    config = joinpath(OUTPUT, "config.toml")

    if !isfile(config)
        @error "$config not found."
        exit(1)
    end

    configfile = TOML.parsefile(config)
    LFP = configfile["lfp"]["save_lfp"] == "n"
    fmax = configfile["lfp"]["fmax"]
    rate = configfile["lfp"]["rate"]

    detect = configfile["detection"]["detect"] == "y"
    fmin_d = configfile["detection"]["fmin"]
    fmax_d = configfile["detection"]["fmax"]
    stdmin = configfile["detection"]["stdmin"]
    stdmax = configfile["detection"]["stdmax"]
    event = configfile["detection"]["type"] == "both" ? 0 : configfile["detection"]["type"] == "pos" ? 1 : -1
    ref = configfile["detection"]["refractory"]          # ms

    shapes = configfile["sorting"]["save_shapes"] == "n"
    fmin_s = configfile["sorting"]["fmin"]
    fmax_s = configfile["sorting"]["fmax"]
    dpre = configfile["sorting"]["delta_pre"]             # ms
    dpost = configfile["sorting"]["delta_post"]            # ms
    factor = configfile["sorting"]["interpolation_factor"]
    spline = configfile["sorting"]["cubic_spline"] == "y"

    #----------------------------------------------

    meta = joinpath(OUTPUT, "meta.toml")
    if !isfile(meta)
        @error "$meta not found. Something must have gone wrong with the bash script."
        exit(1)
    end

    infofile = TOML.parsefile(meta)
    type = infofile["info"]["type"]          # Type of data
    Nchans = infofile["info"]["Nchans"]     # Number of channels
    Nsamples = infofile["info"]["Nsamples"]   # Number of samples per channel

    if type == "MCS"
        
        Tick = infofile["chans"]["Tick"]      # Sampling interval in us
        ADZero = infofile["chans"]["ADZero"]    # ADZero
        ConversionFactor = infofile["chans"]["ConversionFactor"] # D/A conversion factor
        Exponent = infofile["chans"]["Exponent"]  # Exponent for the conversion factor
        #----------------------------------------------
        srate = 1E6 / Tick                          # Sampling rate in Hz
        c = ConversionFactor * 10.0^Exponent  # Conversion from AD to physical units (V)
        d = -ADZero * c                       # Conversion from AD to physical units (V)
        #----------------------------------------------
    elseif type == "3BRAIN"
        srate = infofile["info"]["FrameRate"]          # Sampling rate in Hz
       
        # These could be excluded, 3BRAINS does not have these values
        # but the function still wants to output them
        #----------------------------------------------
        Tick = 0
        ADZero = 0
        ConversionFactor = 0 # D/A conversion factor
        Exponent = 0
        #----------------------------------------------

        maxA = infofile["info"]["MaxAnalogValue"]
        minA = infofile["info"]["MinAnalogValue"]
        maxD = infofile["info"]["MaxDigitalValue"]
        minD = infofile["info"]["MinDigitalValue"]
        c = (maxA - minA) / (maxD - minD)
        d = minA
    end
    # Convert fmin_d, fmax_d, fmin_s, fmax_s, dpre, dpost to Float32, before calling the functions
    bpfilt = SpQEphysTools.prepare_bandpass(Float32(fmin_d), Float32(fmax_d), Float32(srate))
    bpfilt_s = SpQEphysTools.prepare_bandpass(Float32(fmin_s), Float32(fmax_s), Float32(srate))
    lpfilt = SpQEphysTools.prepare_lowpass(Float32(fmax), Float32(srate))

    #-- Settings struct ---------------------------
    return SpQEphysTools.settings(OUTPUT, LFP, type, fmax, rate, detect, fmin_d, fmax_d, stdmin, stdmax, event, ref, shapes, fmin_s, fmax_s, dpre, dpost, factor, spline, Nchans, Nsamples, Tick, ADZero, ConversionFactor, Exponent, srate, c, d, lpfilt, bpfilt, bpfilt_s)

end # function parse_toml_files ----------------


