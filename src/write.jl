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

function tideup_output(OUTPUT)
    # If s.LFP is true, we tide up the OUTPUT folder by moving
    # the LFP files (lfp_*.jld2) to a subfolder /LFP
    if s.LFP
        mkpath(OUTPUT * "/LFP")
        #run(`mkdir -p $OUTPUT/LFP`)
        for file in readdir(OUTPUT) # Loop over the files in the OUTPUT folder
            if startswith(file, "lfp_") && endswith(file, ".jld2")
                #run(`mv $OUTPUT/$file $OUTPUT/LFP/`)
                mv(OUTPUT * "/" * file, OUTPUT * "/LFP/" * file)
            end
        end
    end

    # If s.detect is true, we tide up the OUTPUT folder by moving
    # the SPK files (spk_*.txt) to a subfolder /SPK
    if s.detect
        mkpath(OUTPUT * "/SPK")
        #run(`mkdir -p $OUTPUT/SPK`)
        for file in readdir(OUTPUT) # Loop over the files in the OUTPUT folder
            if startswith(file, "spk_") && endswith(file, ".txt")
                #run(`mv $OUTPUT/$file $OUTPUT/SPK/`)
                mv(OUTPUT * "/" * file, OUTPUT * "/SPK/" * file)
            end
        end
    end

    # If s.shapes is true, we tide up the OUTPUT folder by moving
    # the WAV files (wav_*.jld2) to a subfolder /WAV
    if s.shapes
        #run(`mkdir -p $OUTPUT/WAV`)
        mkpath(OUTPUT * "/WAV")
        for file in readdir(OUTPUT) # Loop over the files in the OUTPUT folder
            if startswith(file, "wav_") && endswith(file, ".jld2")
                #run(`mv $OUTPUT/$file $OUTPUT/WAV/`)
                mv(OUTPUT * "/" * file, OUTPUT * "/WAV/" * file)
            end
        end
    end


    if s.detect
        all_spk = []
        tmq = true
        for i in 1:SpQEphysTools.s.Nchans
            if isfile(OUTPUT * "/SPK/spk_$i.txt") && filesize(OUTPUT * "/SPK/spk_$i.txt") > 0
                tmp = readdlm(OUTPUT * "/SPK/spk_$i.txt")
                tmp[:, 2] = tmp[:, 2] .* i
                if tmq
                    all_spk = tmp
                    tmq = false
                else
                    all_spk = vcat(all_spk, tmp)
                end
            end
        end
        all_spk = sortslices(all_spk, dims=1)
        writedlm(OUTPUT * "/spk.txt", all_spk)
    end

end
