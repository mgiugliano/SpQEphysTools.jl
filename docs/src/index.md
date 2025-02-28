```@meta
CurrentModule = SpiQ
```

# SpiQ

Documentation for [SpiQ](https://github.com/mgiugliano/SpiQ.jl).

List of functions:
- extract_peaks(xf::Array{Float32,1}, thr::Float32, dpre::Float32, dpost::Float32, ref::Float32, event::Int64, srate::Float32)
- prepare_bandpass(lowcut::Float32, highcut::Float32, fs::Float32)::ZPG
- bandpass(data::Array{Float32,1}, filt::ZPG)
- prepare_lowpass(cutoff::Float32, fs::Float32)
- lowpass_and_dec(data::Array{Float32,1}, filt::ZPG, rate::Float32, fs::Float32)
- allocate_Float32vector(num_elements::Int64)
- meminfo_julia()
- preproc_chan(fname::String, chan::Int, datasetname::String)
- parse_toml_files(OUTPUT)::settings
- tideup_output(OUTPUT)
- report(text::String, filename::String)
- active_el(filename::String, threshold::Float64)::Int
- extract_bursts(pathname::String, nActiveEl::Int, ref::Float64)::Int
- plot_raster(pathname::String, fraction::Float64=1.0, title::Bool=false)::Int
- plot_frequencies(pathname::String, fraction::Float64=1.0, title::Bool=false)::Int



```@index
```

```@autodocs
Modules = [SpiQ]
```
