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

module SpQEphysTools

# -- import statement(s) --
using TOML
using HDF5
using DSP
using Statistics
using JLD2
using DelimitedFiles
using Parameters
using .Sys

# -- type definition(s) --
include("types.jl")

# -- global variable(s) --
s = settings()

# -- include statement(s) --
include("memory.jl")
include("readTOML.jl")
include("filtering.jl")
include("detect.jl")
include("write.jl")
include("process_chan.jl")

end # module
