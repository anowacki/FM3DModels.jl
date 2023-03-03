# FM3DModels.jl

FM3DModels.jl is a Julia package to read models and some other file
formats used and written by Nick Rawlinson's FM3D and FMTOMO software.

This is beta software and there may be bugs or incorrect documentation.
I have published this in case it is useful, but I do not guarantee to
continue to maintain this software.

Bug reports and pull requests for new functionality are welcome.

## Installation

```julia
julia> import Pkg

julia> import Pkg.add(url="https://github.com/anowacki/FM3DModels.jl")
```

## Usage

Docstrings describe the purpose and use

FM3DModels exports the `FM3DModel` type, which holds all the information
about an FM3D model.

### Reading a model
```julia
julia> using FM3DModels

julia> mod = FM3DModel("path/to/vgrids.in"); # Read a FM3D grid from disk
```

### Evaluating a model at a given point
```julia
julia> FM3DModels.interpolate(model, lon, lat, radius)
```
