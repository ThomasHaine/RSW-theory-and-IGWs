using JLD2, FileIO
using Printf
using Makie, GLMakie
using Interpolations
using BenchmarkTools

# This line fails because package Oceananigans is needed, but it conflicts with other package in this environment
filename = "../../Oceananigans/RSW_channel_adjustment/RSW_geostrophic_adjustment.jld2"
#model_data = load(filename)

# Therefore load individual parts of the Oceananigans datafile
timestep = 1300
timestep = 0
fld_key = "timeseries/η/"*string(timestep)
time_key = "timeseries/t/"*string(timestep)
η = load(filename,fld_key)
η = η[:,:,1]
t = load(filename,time_key)

function load_grid(type,supscript)
grid = load(filename,"grid/"*type*supscript)
grid_L = load(filename,"grid/L"*type)
grid_N = load(filename,"grid/N"*type)
grid_N_halo = Int((size(grid,1) - grid_N)/2)
@printf("size(grid) = [%d], grid_N = [%d], halo = [%d]\n",size(grid,1),grid_N,grid_N_halo)
grid = grid[grid_N_halo+1:size(grid,1)-grid_N_halo]
return grid, grid_L
end

# Grid coordinates at cell centers:
grid_x, Lx = load_grid("x","ᶜᵃᵃ")
grid_y, Ly = load_grid("y","ᵃᶜᵃ")

# plots
kwargs = (
          fillrange = true,
        levels = range(-2.2,2.2,length=64),
      colormap = :seismic
    )
fig,ax,plt = contourf(grid_x / Lx, grid_y / Ly, η, figure = (resolution = (700, 600),),axis = (title = "Test", xlabel = L"x/L_x",ylabel = L"y/L_x",); kwargs...)
ax.aspect = DataAspect()
Colorbar(fig[1,2],plt)
display(fig)

interp_linear = LinearInterpolation((grid_x/Lx, grid_y/Ly), η, extrapolation_bc = Flat())

@benchmark interp_linear(0.1,-0.3)