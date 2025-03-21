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
    extract_peaks(xf::Array{Float32,1}, thr::Float32, dpre::Float32, dpost::Float32, ref::Float32, event::Int64, srate::Float32)

    Performs a threshold crossing detection on a given signal (xf), and extracts the peaks of the signal in terms of their
    index/time of occurrence and "polarity" (i.e., sign of the peak). The function is designed to work with a single channel
    signal, and it is used in the context of spike detection in electrophysiological recordings.

### Arguments
    It requires the following parameters:
    - `xf::Array{Float32,1}` : the signal to be analyzed.
    - `thr::Float32` : the threshold value for the detection.
    - `dpre::Float32` : the duration of the pre-event window (in ms), for storing the event waveform.
    - `dpost::Float32` : the duration of the post-event window (in ms), for storing the event waveform.
    - `ref::Float32` : the refractory period (in ms), to avoid multiple detections of the same event.
    - `event::Int64` : the type of event to be detected (1 for positive, -1 for negative, 0 for both).
    - `srate::Float32` : the sampling rate of the signal (in Hz).

### Author(s)
- Michele Giugliano - michele.giugliano@unimore.it
"""
function extract_peaks(xf::Array{Float32,1}, thr::Float32, dpre::Float32, dpost::Float32, ref::Float32, event::Int64, srate::Float32)
    idx = Vector{Vector{Float32}}()             # Index of events and polarity

    wpre = Int64(ceil(1e-3 * dpre * srate))   # Pre window, in samples
    wpost = Int64(ceil(1e-3 * dpost * srate))   # Post window, in samples
    rref = Int64(ceil(1e-3 * ref * srate))    # Refr/dead period, in samples
    href = Int64(ceil(1e-3 * 0.5 * ref * srate)) # Half refr/dead period, in samples

    if event == 1       # "positive" threshold crossings
        y = findall(xf[wpre+2:end-wpost-2] .> thr) .+ (wpre + 1) # indx of elms > thr
    elseif event == -1  # "negative" threshold crossings
        y = findall(xf[wpre+2:end-wpost-2] .< -thr) .+ (wpre + 1) # indx of elms < thr
    else                # "both" threshold crossings
        y = findall(abs.(xf[wpre+2:end-wpost-2]) .> thr) .+ (wpre + 1) # indx of elms > thr
    end

    last = 0                          # index last event detected so far
    for i in eachindex(y)               # over all threshold crossings (indexes)
        if y[i] >= last + rref         # current event, after last refractory/dead period
            A = xf[y[i]:y[i]+href-1]   # extract n samples = to half refr/dead period (volt)
            A = abs.(A)                # take their abs (so I can always use "maximum" below)
            iaux = findall(A .== maximum(A)) # introduces alignment: takes indx of max (index)
            index = y[i] + iaux[1] - 1  # index of max value in original signal (xf)
            polar = sign(xf[index])     # polarity of max value (xf)
            #append!(idx, iaux .+ (y[i] -1)); # append index of max to events list (index)
            push!(idx, [index, polar])  # append index of max to events list (index)
            last = index               # update index of last event detected so far
        end
    end

    if isempty(idx)                    # if no events detected
        return []                      # return empty list
    else
        return reduce(vcat, idx')      # return the list of events (index and polarity)
    end
end # extract_peaks ----------------------------



