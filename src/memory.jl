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

"""
    allocate_Float32vector(num_elements::Int64)

Allocate a Float32 vector of a given size, handling OutOfMemoryError exceptions
"""
function allocate_Float32vector(num_elements::Int64)
  try
    # Attempt to allocate the vector
    vec = Vector{Float32}(undef, num_elements)
    #@info "Memory successfully allocated ($(sizeof(vec)/1e6) MB)."
    return vec
  catch e
    if isa(e, OutOfMemoryError)
      @error "ERROR: Not enough memory available."
      # Handle the error appropriately (e.g., exit, try smaller size)
      exit(1)
    else
      rethrow(e)  # Re-throw other exceptions
    end
  end
end # allocate_Float32vector ----------------------



function meminfo_julia()
  # @printf "GC total:  %9.3f MiB\n" Base.gc_total_bytes(Base.gc_num())/2^20
  # Total bytes (above) usually underreports, thus I suggest using live bytes (below)
  #@info "GC live:   %9.3f MiB\n" Base.gc_live_bytes()/2^20
  #@info "JIT:       %9.3f MiB\n" Base.jit_total_bytes()/2^20
  #@info "Max. RSS:  %9.3f MiB\n" Sys.maxrss()/2^20
  println(Sys.maxrss()/2^30)
  end

