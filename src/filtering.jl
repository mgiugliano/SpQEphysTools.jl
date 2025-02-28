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

"""
    prepare_bandpass(lowcut::Float32, highcut::Float32, fs::Float32)

    This function prepares (once for all) the bandpass filter to be used for spike detection/sorting.
    It returns the filter object. It receives the low and high cutoff frequencies (in Hz) and the sampling rate (in Hz).
    Its purpose is to avoid redefining the filter at each call of the bandpass function.
    Its output is used as input to the bandpass function.
"""
function prepare_bandpass(lowcut::Float32, highcut::Float32, fs::Float32)::ZPG
    nyquist = 0.5 * fs                 # Nyquist frequency
    low = lowcut / nyquist             # Normalized cutoff frequencies
    high = highcut / nyquist           # Normalized cutoff frequencies
    TYPE = Bandpass(low, high)         # Bandpass filter
    DESIGN = Elliptic(2, 0.1, 40)      # 2nd order, 0.1 dB passband ripple, 40 dB stopband attenuation
    return digitalfilter(TYPE, DESIGN) # Return the filter object (ZPG - ZeroPoleGain)
end # prepare_bandpass ----------------------------

"""
    bandpass(data::Array{Float32,1}, lowcut::Float32, highcut::Float32, fs::Float32)

    This function applies a bandpass filter to the input data. It receives the data (as a 1D array of Float32),
    and the filter parameters previously defined by prepare_bandpass. It returns the filtered data.
"""
@inline function bandpass(data::Array{Float32,1}, filt::ZPG)
    return filtfilt(filt, data)
end # bandpass -----------------------------------


"""
    prepare_lowpass(cutoff::Float32, fs::Float32)

    This function prepares (once for all) the lowpass filter to be used for LFP extraction.
    It returns the filter object. It receives the low cutoff frequency (in Hz) and the sampling rate (in Hz).
    Its purpose is to avoid redefining the filter at each call of the lowpass function.
    Its output is used as input to the lowpass function.
"""
function prepare_lowpass(cutoff::Float32, fs::Float32)
    nyquist = 0.5 * fs         # Nyquist frequency
    low = cutoff / nyquist     # Normalized cutoff frequency
    TYPE = Lowpass(low)      # Lowpass filter
    DESIGN = Butterworth(5)    # 5th order Butterworth filter
    return digitalfilter(TYPE, DESIGN) # Return the filter object (ZPG - ZeroPoleGain)
end # prepare_lowpass ----------------------------

"""
    lowpass_and_dec(data::Array{Float32,1}, cutoff::Float32, rate::Float32, fs::Float32)

    This function applies a lowpass filter to the input data and further decimates the result.
    It receives the data (as a 1D array of Float32), and the filter parameters previously defined
    by prepare_lowpass. It returns the filtered and decimated data.
"""
#function lowpass_and_dec(data::Array{Float32,1}, cutoff::Float32, rate::Float32, fs::Float32)
@inline function lowpass_and_dec(data::Array{Float32,1}, filt::ZPG, rate::Float32, fs::Float32)
    xf = filtfilt(filt, data)       # Apply lowpass filter
    decimate = Int(ceil(fs / rate)) # Decimation factor
    return xf[1:decimate:end]       # Return the filtered and decimated data
end # lowpass_and_dec ----------------------------


