"""
### FM3DModels

Read the `vgrid.in`-type files used by Nick Rawlinson's fast marching code
`FM3D` and used by his tomography code `FMTOMO`.
"""
module FM3DModels

using DelimitedFiles: readdlm
using Statistics: mean

using Interpolations: AbstractInterpolation, CubicSplineInterpolation, Flat
using WriteVTK: vtk_grid, vtk_point_data


export
    FM3DModel,
    layer_perturbation!,
    layer_perturbation,
    modify_grid,
    reference_perturbation!,
    reference_perturbation,
    trim,
    write_vtk

"""
    FM3DModel

Struct holding a velocity grid in geographic coordinates.
`dlon`, `dlat` and `dr` are grid spacings in °, ° and km, respectively.
`lon0`, `lat0` and `r0` are the grid origin (°, ° and km), and
`nlon`, `nlat` and `nr` are the number of grid points in each
dimension.

The vectors `lon`, `lat` and `r` hold the coordinates of the grid points
in each dimension, whilst `grid` hold the values of the grid,
accessed by `grid[ilon,ilat,ir]`; i.e., longitude varies in the first
dimension, latitude in the second and radius in the last.

!!! note
    The `grid` contained within an `FM3DModel` is mutable, but changing
    `grid` is not supported as the interpolation is no longer valid.
    Instead, use [`modify_grid`](@ref).
"""
struct FM3DModel{T, R<:AbstractRange, I<:AbstractInterpolation}
    nlon::Int
    nlat::Int
    nr::Int
    dlon::T
    dlat::T
    dr::T
    lon0::T
    lat0::T
    r0::T
    lon::R
    lat::R
    r::R
    grid::Array{T,3}
    itp::I
end

"""
    FM3DModel(file) -> ::FM3DModel

Construct a `FM3DModel` by reading from a `file` in the `vgrids.in` format.
"""
FM3DModel(file) = read_fm_model(file)

"""
    FM3DModel(lon, lat, r, grid)

Construct a `FM3DModel` with a `lon`gitudes, `lat`itudes (°) and
`r`adii (km).  `grid` is the velocity grid which should have the same number
of points in each dimension as the coordinates vectors/ranges.

Computation of the grid spacing and number of points is done simply
by looking at the length of the coordinate vectors and the spacing
of the first two coordinates, if they are vectors, or their step if
they are ranges.
"""
function FM3DModel(lon::AbstractRange, lat::AbstractRange, r::AbstractRange, grid)
    T = promote_type(eltype.((lon, lat, r, grid))...)
    nlon, nlat, nr = length.((lon, lat, r))
    lon0, lat0, r0 = T.(first.((lon, lat, r)))
    dlon, dlat, dr = T.(step.((lon, lat, r)))
    itp = CubicSplineInterpolation((lon, lat, r), grid, extrapolation_bc=Flat())
    FM3DModel(nlon, nlat, nr, dlon, dlat, dr, lon0, lat0, r0,
        lon, lat, r, grid, itp)
end

_step(x::AbstractRange) = step(x)
_step(x) = x[2] - x[1]

"""
    grid_index(m, ilon, ilat, ir) -> i

Return the global grid index for a velocity point within the model `m`
with indices in longitude, latitude and radius respectively of `ilon`,
`ilat` and `ir`.
"""
grid_index(m::FM3DModel, ilon, ilat, ir) =
    LinearIndices((m.nr, m.nlat, m.nlon))[ir,ilat,ilon]

"""
    grid_indices(m, i) -> ilon, ilat, ir

Return the coordinate indices for the global grid index `i`.
"""
grid_indices(m::FM3DModel, i) = reverse(Tuple(CartesianIndices((m.nr, m.nlat, m.nlon))[i]))

"""
    read_fm_model(file) -> x, y, z, grid

Read an FM3D model file.  This has the format of a FM3D 'vgrids.in'-type file, i.e.:

```
1 1            # Number of grids and number of velocity types
11 21 21       # Number grid points in radius, latitude, longitude (nr, nlat, nlon)
10 0.001 0.001 # Radial grid spacing (km) and lat, long spacing (radians) (dr, dlat, dlon)
6270 0.78 0.78 # Grid origin in radius (km), lat and long (radians) (r0, lat0, lon0)
grid[1,1,1]
grid[2,1,1]
⋮
grid[nlon-1,nlat,nr]
grid[nlon,nlat,nr]
```

`grid` in the input file is defined to loop fastest over longitude and
slowest over radius.  As returned, `grid` has coordinates `(lon, lat, r)`.
"""
function read_fm_model(file)
    d = readdlm(file)
    d[1,1:2] == [1, 1] ||
        warn("File '$file' contains more than one grid; only reading the first")
    nz, ny, nx = Int.(d[2,1:3])
    dz, dy, dx = Float64.(d[3,1:3])
    z0, y0, x0 = Float64.(d[4,1:3])
    x0, y0, dx, dy = rad2deg.((x0, y0, dx, dy))
    x = range(x0, step=dx, length=nx)
    y = range(y0, step=dy, length=ny)
    z = range(z0, step=dz, length=nz)
    grid = reshape(Float64.(d[5:end,1]), nx, ny, nz)
    FM3DModel(x, y, z, grid)
end

"""
    interpolate(m::FM3DModel; n, d) -> m′

Return a copy of `m` where the grid spacing is reduced by `n` times, or
where the new grid spacing is `d` km in each direction.
"""
function interpolate(m::FM3DModel{T}; n=nothing, d=nothing)  where T
    if n != nothing
        n1, n2, n3 = length(n) == 3 ? n : (n, n, n)
        nlon, nlat, nr = (n1*(m.nlon - 1) + 1), (n2*(m.nlat - 1) + 1), (n3*(m.nr - 1) + 1)
        lon = range(m.lon0, stop=m.lon[end], length=nlon)
        lat = range(m.lat0, stop=m.lat[end], length=nlat)
        r = range(m.r0, stop=m.r[end], length=nr)
    elseif d != nothing
        (d1, d2, d3) = length(d) == 3 ? d : (d, d, d)
        lon = range(m.lon0, stop=m.lon[end], step=d1)
        lat = range(m.lat0, stop=m.lat[end], step=d2)
        r = range(m.r0, stop=m.r[end], step=d3)
    else
        throw(ArgumentError("must specify one of `d` or `n`"))
    end
    nlon, nlat, nr = length.((lon, lat, r))
    grid = Array{T}(undef, nlon, nlat, nr)
    m′ = FM3DModel(lon, lat, r, grid)
    for (ir, rr) in enumerate(r), (ilat, la) in enumerate(lat), (ilon, lo) in enumerate(lon)
        grid[ilon,ilat,ir] = interpolate(m, lo, la, rr)
    end
    m′
end

"""
    interpolate(m::FM3DModel, lon, lat, r) -> val

Return the value of the model `m` interpolated to point `lon`, `lat`, `r`,
calculated using cubic spline interpolation.
"""
interpolate(m::FM3DModel, lon, lat, r) = m.itp(lon, lat, r)

"""
    trim(m::FM3DModel; lon=(-Inf,Inf), lat=(-Inf,Inf), r=(-Inf,Inf), R=6371, dep=R.-r) -> m′

Cut out a section of model `m` where points lie within the `lon`gitude,
`lat`itude (both °) and `r`adius bounds (km) provided as keyword
arguments.  If `dep` is provided instead of `r`, then trim based on a depth
range (km), where `R` specifies the Earth radius in km.
"""
function trim(m::FM3DModel;
        lon=(-Inf,Inf), lat=(-Inf,Inf), r=(-Inf,Inf), R=6371, dep=nothing)
    ilon = lon[1] .<= m.lon .<= lon[end]
    ilat = lat[1] .<= m.lat .<= lat[end]
    if dep !== nothing
        r = R .- reverse(dep)
    end
    ir = r[1] .<= m.r .<= r[end]
    # Get ranges back
    lons = m.lon[findfirst(ilon):findlast(ilon)]
    lats = m.lat[findfirst(ilat):findlast(ilat)]
    rs = m.r[findfirst(ir):findlast(ir)]
    FM3DModel(lons, lats, rs, m.grid[ilon,ilat,ir])
end

"""
    layer_perturbation(m::FM3DModel) -> m′

Return a copy of `m` where each radial layer in `m.grid`
is replaced with the percentage perturbation from the layer mean.
"""
function layer_perturbation!(m::FM3DModel)
    modify_grid(m) do grid
        grid′ = similar(grid)
        for i in axes(grid, 3)
            layer = @view(grid[:,:,i])
            layer_mean = mean(layer)
            grid′[:,:,i] .= 100*(layer .- layer_mean)./layer_mean
        end
        grid′
    end
end

"""
    reference_perturbation(m, ref) -> m′

Convert model `m` into a relative model with respect to the reference model
`ref`.  Both models must have the exact same grid layout.
"""
function reference_perturbation(m::FM3DModel, ref::FM3DModel)
    for field in (:lon, :lat, :r)
        getfield(m, field) ≈ getfield(ref, field) ||
            throw(ArgumentError("model and reference model have different layouts"))
    end
    modify_grid(m) do grid
        100 .* (grid .- ref.grid)./ref.grid
    end
end

"""
    modify_grid(f, m) -> m′

Return a copy of the `FM3DModel` `m` with its grid altered by the
function `f`, which takes the grid as its single argument and
should return a grid of the correct dimensions.

!!! note
    The grid of `m` is passed to `f`, so make sure to copy it
    first if you plan on modifying it in place before returning in `f`.

# Example

Create a model then add one to each grid point, creating a new model.

```
julia> m = FM3DModel(1:3, 1:3, -2:0, zeros(3, 3, 3));

julia> extrema(m.grid)
(0.0, 0.0)

julia> m′ = modify_grid(m) do grid
           grid .+ 1
       end;

julia> extrema(m′.grid)
(1.0, 1.0)
```

In this case, one could also simply call `m′ = modify_grid(x->x.+1, m)`.
"""
function modify_grid(f, m::FM3DModel)
    grid = f(m.grid)
    FM3DModel(m.lon, m.lat, m.r, grid)
end

"""
    write_vtk(m::FM3DModel, filebase, name="Velocity";
        lon=(-Inf,Inf), lat=(-Inf,Inf), r=(-Inf,Inf), r0=0, layer_mean=false) -> filename

Write a `FM3DModel` as a VTK file called `filename`, given the file name
stem `filebase`.  Optionally specify the name of the velocity field as `name`.

Limit the saved region by supplying a tuple for one or mor of `lon`, `lat`
and `r`.  To shift the radial origin to be other than that in the model
(for instance to plot relative to sea level) set `r0`, and radii will have
this value subtracted.

If `layer_mean` is `true`, then subtract the mean value from each radial layer.

Longitude and latitude are converted to km by using the middle radial value
and assuming a spherical geometry.  The mean latitude is used to scale
longitude to km at this radius.  Note therefore that this will introduce
significant distortion at high latitudes.
"""
function write_vtk(m::FM3DModel, filebase, name="Velocity";
        lon=(-Inf,Inf), lat=(-Inf,Inf), r=(-Inf,Inf), r0=0,
        layer_mean=false)
    layer_mean && (m = layer_perturbation(m))
    deg2km_lat = dist_per_degree_arc(mean(m.r))
    deg2km_lon = deg2km_lat*cosd(mean(m.lat))
    ilon = lon[1] .<= m.lon .<= lon[end]
    ilat = lat[1] .<= m.lat .<= lat[end]
    ir = r[1] .<= m.r .<= r[end]
    x = deg2km_lon.*(m.lon[ilon] .- m.lon0)
    y = deg2km_lat.*(m.lat[ilat] .- m.lat0)
    z = m.r[ir] .- r0
    filename = vtk_grid(filebase, x, y, z) do vtk
        vtk[name] = m.grid[ilon,ilat,ir]
    end
    if length(filename) > 1
        @warn("More than one output file written by `vtk_grid`")
        filename
    else
        filename[1]
    end
end

"""
    read_frechet_derivs(m, file)

Read a set of Fréchet derivatives for an `FM3DModel` `m`.  The `file`
(usually called `frechet.dat`) is assumed to match the dimensions of
model `m`.

!!! note
    At present, only the sum of the absolute values of the Fréchet derivatives
    at each grid point are read in.  There is no information about which rays
    they correspond to.
"""
function read_frechet_derivs(m::FM3DModel{T}, file; normalise=false) where T
    grid = zeros(T, axes(m.grid))
    open(file, "r") do f
        while !eof(f)
            tokens = _get_tokens(readline(f))
            length(tokens) != 2 && continue
            i = parse(Int, tokens[1])
            dv = parse(T, tokens[2])
            ilon, ilat, ir = grid_indices(m, i)
            grid[ilon,ilat,ir] += abs(dv)
        end
    end
    normalise && (grid ./= maximum(grid))
    FM3DModel(m.lon, m.lat, m.r, grid)
end

"""
    read_rays(file) -> lons, lats, rs

Read the ray paths used by FM3D from `file`, usually called `rays.in`,
returning vectors of vectors of `lon`gitude, `lat`itude (both °) and `r`adius (km).
"""
function read_rays(file)
    lons = Vector{Float64}[]
    lats = similar(lons)
    rs = similar(lons)
    iline = 1
    nrays = 0
    num_ray_sections = 0
    iray_section = 0

    open(file, "r") do f
        expected_line = :header
        while !eof(f)
            line = readline(f)
            tokens = _get_tokens(line)
            if expected_line === :header
                length(tokens) == 5 || _line_read_error(iline, line)
                receiver_id, source_id, path_id = parse.(Int, tokens[1:3])
                is_normal_path = parse(Int, tokens[4]) == 0
                num_ray_sections = parse(Int, tokens[5])
                if num_ray_sections == 0
                    expected_line = :header
                    iline += 1
                else
                    iray_section = 1
                    push!(lons, [])
                    push!(lats, [])
                    push!(rs, [])
                    expected_line = :ray_header
                    iline += 1
                end
            elseif expected_line === :ray_header
                length(tokens) == 4 || _line_read_error(iline, line)
                num_points, region_id = parse.(Int, tokens[1:2])
                has_diffraction, has_headwave = tokens[3:4] .== "F"
                for _ in 1:num_points
                    line = readline(f)
                    tokens = _get_tokens(line)
                    length(tokens) == 3 || _line_read_error(iline, line)
                    _r, _lat, _lon = parse.(Float64, tokens)
                    push!(lons[end], rad2deg(_lon))
                    push!(lats[end], rad2deg(_lat))
                    push!(rs[end], _r)
                end
                iline += num_points + 1
                if iray_section == num_ray_sections
                    expected_line = :header
                else
                    expected_line = :ray_header
                end
            end
        end
    end
    lons, lats, rs
end

"""
    read_receivers(file) -> lon, lat, dep

Read the receivers used by FM3D from `file`, usually called `receivers.in`,
returning vectors of `lon`gitude and `lat`itude (°) and `dep`th (km).
"""
function read_receivers(file)
    lons = Float64[]
    lats = similar(lons)
    deps = similar(lons)
    iline = 1
    nreceivers = 0

    open(file, "r") do f
        expected_line = :header
        while !eof(f)
            line = readline(f)
            tokens = _get_tokens(line)
            if expected_line === :header
                iline == 1 || error("header line on wrong line ($line)")
                length(tokens) == 1 || _line_read_error(iline, line)
                nreceivers = parse(Int, tokens[1])
                expected_line = :coordinates
                iline += 1
            elseif expected_line === :coordinates
                length(tokens) == 3 || _line_read_error(iline, line)
                dep, lat, lon = parse.(Float64, tokens[1:3])
                push!(lons, lon)
                push!(lats, lat)
                push!(deps, dep)
                expected_line = :number_of_paths
                iline += 1
            elseif expected_line === :number_of_paths
                length(tokens) == 1 || _line_read_error(iline, line)
                npaths = parse(Int, tokens[1])
                nlines_to_skip = 2*npaths
                for _ in 1:nlines_to_skip
                    readline(f)
                end
                expected_line = :coordinates
                iline += nlines_to_skip + 1
            end
        end
    end
    if length(lons) != nreceivers
        @warn("Read $(length(lons)) receivers but header specified $nreceivers")
    end
    lons, lats, deps
end

"""
    read_sources(file) -> lon, lat, dep

Read the sources used by FM3D frm `file`, usually called `sources.in`,
returning vectors of `lon`gitude and `lat`itude (°) and `dep`th (km).
"""
function read_sources(file)
    lons = Float64[]
    lats = similar(lons)
    deps = similar(lons)
    iline = 1
    nsources = 0

    open(file, "r") do f
        expected_line = :header
        while !eof(f)
            line = readline(f)
            tokens = _get_tokens(line)
            if expected_line === :header
                iline == 1 || error("header line on wrong line ($iline)")
                length(tokens) == 1 || _line_read_error(iline, line)
                nsources = parse(Int, tokens[1])
                expected_line = :teleseismic
                iline += 1
            elseif expected_line === :teleseismic
                length(tokens) == 1 || _line_read_error(iline, line)
                isteleseismic = parse(Int, tokens[1]) == 1
                if isteleseismic
                    phase = strip(readline(f))
                    iline += 1
                end
                expected_line = :coordinates
                iline += 1
            elseif expected_line === :coordinates
                length(tokens) in (3, 4) || _line_read_error(iline, line)
                dep, lat, lon = parse.(Float64, tokens[1:3])
                push!(lons, lon)
                push!(lats, lat)
                push!(deps, dep)
                expected_line = :number_of_paths
                iline += 1
            elseif expected_line === :number_of_paths
                length(tokens) == 1 || _line_read_error(iline, line)
                npaths = parse(Int, tokens[1])
                nlines_to_skip = 3*npaths
                for _ in 1:nlines_to_skip
                    readline(f)
                end
                expected_line = :teleseismic
                iline += nlines_to_skip + 1
            end
        end
    end
    if length(lons) != nsources
        @warn("Read $(length(lons)) sources but header specified $nsources")
    end
    lons, lats, deps
end

_get_tokens(line) = split(replace(strip(line), r"\s+"=>" "))
_line_read_error(iline, line) = error("Wrong number of tokens on line $iline: \"$(line)\"")

"""Return the distance swept by a degree of arc around a circle at radius `r`."""
dist_per_degree_arc(r) = 2π*r/360

end # module

