# SpQEphysTools - QSpike Tools reinvented - electrophysiology extracellular multichannel batch and parallel preprocessor
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

"""
    extract_peaks(xf::Array{Float32,1}, thr::Float32, dpre::Float32, dpost::Float32, ref::Float32, event::Int64, srate::Float32)

TBW
"""
function extract_peaks(xf::Array{Float32,1}, thr::Float32, dpre::Float32, dpost::Float32, ref::Float32, event::Int64, srate::Float32)
        idx = Vector{Vector{Float32}}()             # Index of events and polarity

        wpre  = Int64(ceil(1e-3 * dpre  * srate));   # Pre window, in samples
        wpost = Int64(ceil(1e-3 * dpost * srate));   # Post window, in samples
        rref  = Int64(ceil(1e-3 * ref  * srate));    # Refr/dead period, in samples
        href  = Int64(ceil(1e-3 * 0.5 * ref  * srate)); # Half refr/dead period, in samples

        if event == 1       # "positive" threshold crossings
            y = findall(xf[wpre+2:end-wpost-2]       .> thr) .+ (wpre+1) # indx of elms > thr
        elseif event == -1  # "negative" threshold crossings
            y = findall(xf[wpre+2:end-wpost-2]      .< -thr) .+ (wpre+1) # indx of elms < thr
        else                # "both" threshold crossings
            y = findall(abs.(xf[wpre+2:end-wpost-2]) .> thr) .+ (wpre+1) # indx of elms > thr
        end

	    last = 0;                          # index last event detected so far
        for i in eachindex(y)               # over all threshold crossings (indexes)
	        if y[i] >= last + rref         # current event, after last refractory/dead period
                A = xf[y[i]:y[i]+href-1]   # extract n samples = to half refr/dead period (volt)
                A = abs.(A)                # take their abs (so I can always use "maximum" below)
                iaux = findall(A .== maximum(A)) # introduces alignment: takes indx of max (index)
                index = y[i] + iaux[1] - 1  # index of max value in original signal (xf)
                polar = sign(xf[index])     # polarity of max value (xf)
                #append!(idx, iaux .+ (y[i] -1)); # append index of max to events list (index)
                push!(idx, [index, polar])  # append index of max to events list (index)
                last = index;               # update index of last event detected so far
	        end
	    end

        return reduce(vcat, idx')           # return the list of events (index and polarity)
    end # extract_peaks ----------------------------



