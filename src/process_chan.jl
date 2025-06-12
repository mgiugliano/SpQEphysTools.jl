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

    c = Float32(s.c)                                            # Conversion factor (D/A)
    d = Float32(s.d)                                            # Conversion factor (D/A)
    if s.type == "MCS"
        REQ = size(data, 1) * 32 / 8 / 1e6                                  # Required memory in MB
        #MEM = Sys.free_memory() / 1e6;                                    # Free memory in MB
        @info "Chan $(chan): $(REQ) MB to read."
        DATA = Float32.(c .* data[:, chan] .+ d)                    # This is time consuming
    elseif s.type == "3BRAIN"
        REQ = size(data, 1) * 32 / 8 / 1e6 / s.Nchans                                  # Required memory in MB
        @info "Chan $(chan): $(REQ) MB to read."
        chanDATA = data[chan:4096:end] # Get the data for the channel
        realREQ = size(chanDATA, 1) * 32 / 8 / 1e6                                  # Required memory in MB
        DATA = Float32.(c .* chanDATA .+ d)                    # This is time consuming CHECK CHECK CHECK
    end
   

    # LFP extraction -------------------------------------
    if s.LFP
        if s.type == "3BRAIN"
            xf = allocate_Float32vector(Int64(size(data, 1) / s.Nchans))
        elseif s.type == "MCS"
            xf = allocate_Float32vector(size(data, 1)) 
        end
        lpfilt = s.lpfilt
        xf = lowpass_and_dec(DATA, lpfilt, s.rate, s.srate)
        println("LFP: $(size(xf))")
        # outname = joinpath(s.OUTPUT, "lfp_$(chan).jld2")
        # open(outname, "w") do f
        #     write(f, "lfp", xf)
        # end
        
        outname = joinpath(s.OUTPUT, "lfp_$(chan).txt")
        writedlm(outname, xf) # Adam: temporary solution. JLD2 files were coming out corrupted
      

        @info "LFP: Chan $(chan) done!"
        xf = nothing # Let's free some memory
    end # LFP --------------------------------------------

    # Spike detection ------------------------------------
    if s.detect
        if s.type == "3BRAIN"
            xf = allocate_Float32vector(Int64(size(data, 1) / s.Nchans))
        elseif s.type == "MCS"
            xf = allocate_Float32vector(size(data, 1)) 
        end
       
        bpfilt = s.bpfilt
        xf = bandpass(DATA, bpfilt)

        noise_std = median(abs.(xf) / 0.6745)
        thr = Float32(s.stdmin * noise_std)
        #thrmax = Float32(s.stdmax * noise_std();

        #@info "Chan $(chan) - detecting spikes...";
        idx = extract_peaks(xf, thr, s.dpre, s.dpost, s.ref, s.event, s.srate)
        tspk = idx
        tspk[:, 1] = tspk[:, 1] / s.srate

        outname = joinpath(s.OUTPUT, "spk_$(chan).txt")
        writedlm(outname, tspk)

        noise_outname = joinpath(s.OUTPUT, "noise_$(chan).txt")
        writedlm(noise_outname, noise_std) 
        
        @info "MUA: Chan $(chan) done! $(length(idx)) events detected."
    end # Spike detection --------------------------------



    if s.shapes
        if s.fmin_s != s.fmin_d || s.fmax_s != s.fmax_d
            @warn "Different filters for spike detection and shapes!"
            if s.type == "3BRAIN"
            xf = allocate_Float32vector(Int64(size(data, 1) / s.Nchans))
        elseif s.type == "MCS"
            xf = allocate_Float32vector(size(data, 1)) 
        end
            bpfilt_s = s.bpfilt_s
            xf = bandpass(DATA, bpfilt_s)

            #noise_std = median(abs.(xf) / 0.6745);
            #thr    = stdmin * noise_std;
            #thrmax = stdmax * noise_std;
        end

        wpre = Int32(ceil(1e-3 * s.dpre * s.srate))   # Pre window, in samples
        wpost = Int32(ceil(1e-3 * s.dpost * s.srate))   # Post window, in samples

        outname = joinpath(s.OUTPUT, "wav_shapes_c$(chan).txt")
        #outname = joinpath(s.OUTPUT, "wav_shapes_c$(chan).jld2")
        open(outname, "a") do f
            for i in eachindex(idx[:, 1])
                ap_sample = Int32(ceil(idx[i, 1] * s.srate))
                if ap_sample - wpre > 0 && ap_sample + wpost < length(xf)
                    start_sample = ap_sample - wpre
                    end_sample = ap_sample + wpost
                    wave = xf[start_sample:end_sample]
                    #write(f, "wav_$(i)", wave)
                    writedlm(f, wave',",") # Adam: temporary solution. JLD2 files were coming out corrupted
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


