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
    preproc_chan(fname, chan, datasetname)

TBW
"""
@noinline function preproc_chan(fname::String, chan::Int, datasetname::String)
    file = h5open(fname, "r")                                          # Open .h5 file (on each worker)
    data = file[datasetname]                                           # Get the dataset
    #data = file["Data/Recording_0/AnalogStream/Stream_0/ChannelData"];# Get the dataset
    REQ = size(data, 1) * 32 / 8 / 1e6                                  # Required memory in MB
    #MEM = Sys.free_memory() / 1e6;                                    # Free memory in MB
    @info "Chan $(chan): $(REQ) MB to read."

    c = Float32(s.c)                                            # Conversion factor (D/A)
    d = Float32(s.d)                                            # Conversion factor (D/A)

    DATA = Float32.(c .* data[:, chan] .+ d)                    # This is time consuming

    # LFP extraction -------------------------------------
    if s.LFP
        xf = allocate_Float32vector(size(data, 1))
        lpfilt = s.lpfilt
        xf = lowpass_and_dec(DATA, lpfilt, s.rate, s.srate)
        outname = joinpath(s.OUTPUT, "lfp_$(chan).jld2")
        open(outname, "w") do f
            write(f, "lfp", xf)
        end
        @info "LFP: Chan $(chan) done!"
        xf = nothing # Let's free some memory
    end # LFP --------------------------------------------

    # Spike detection ------------------------------------
    if s.detect
        xf = allocate_Float32vector(size(data, 1))
        bpfilt = s.bpfilt
        xf = bandpass(DATA, bpfilt)

        noise_std = median(abs.(xf) / 0.6745)
        thr = Float32(s.stdmin * noise_std)
        #thrmax = Float32(s.stdmax * noise_std();

        #@info "Chan $(chan) - detecting spikes...";
        idx = extract_peaks(xf, thr, s.dpre, s.dpost, s.ref, s.event, s.srate)

        if length(idx) == 0
            @warn "Chan $(chan) - no spikes detected!"
            return
        end

        tspk = idx
        tspk[:, 1] = tspk[:, 1] / s.srate

        outname = joinpath(s.OUTPUT, "spk_$(chan).txt")
        writedlm(outname, tspk) # Save the spike times, in seconds, to a text file
        @info "MUA: Chan $(chan) done! $(length(idx)) events detected."
    end # Spike detection --------------------------------



    if s.shapes && (length(idx) > 0)
        if s.fmin_s != s.fmin_d || s.fmax_s != s.fmax_d
            @warn "Different filters for spike detection and shapes!"
            xf = allocate_Float32vector(size(data, 1))
            bpfilt_s = s.bpfilt_s
            xf = bandpass(DATA, bpfilt_s)

            #noise_std = median(abs.(xf) / 0.6745);
            #thr    = stdmin * noise_std;
            #thrmax = stdmax * noise_std;
        end

        wpre = Int32(ceil(1e-3 * s.dpre * s.srate))   # Pre window, in samples
        wpost = Int32(ceil(1e-3 * s.dpost * s.srate))   # Post window, in samples

        outname = joinpath(s.OUTPUT, "spk_shapes_c$(chan).jld2")
        open(outname, "w") do f
            for i in eachindex(idx[:, 1])
                if idx[i, 1] - wpre > 0 && idx[i] + wpost < length(xf)
                    wave = xf[idx[i, 1]-wpre:idx[i, 1]+wpost]
                    write(f, "wav_$(i)", wave)
                end # if
            end # for
        end # open
        @info "WAV: Chan $(chan) done!"
    end # Shapes ----------------------------------------

    # Free memory
    DATA = nothing
    data = nothing
    xf = nothing
    tspk = nothing
    #SpQEphysTools.meminfo_julia()

    close(file)     # Close the file
end # preproc_chan
#----------------------------------------------


