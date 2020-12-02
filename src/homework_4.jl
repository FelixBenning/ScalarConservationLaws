### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ ccd42400-3266-11eb-3d2c-97c96caf3e78
begin
	using Printf: @sprintf
	using Zygote
	using Plots: Plots, plot, plot!
	using Roots: find_zero
	using PlutoUI
end

# ╔═╡ 9863f500-3267-11eb-32e7-b5a82fb92dc9
abstract type Junction end

# ╔═╡ 61bc4c92-3268-11eb-1197-7786485b47d1
abstract type Line end

# ╔═╡ e31a36f0-3284-11eb-08a3-77c7fc3c47e6
abstract type Network end

# ╔═╡ 663ebce2-3267-11eb-1394-318abf53bd13
mutable struct ProductionLine <: Line
	flux
	capacity
	velocity
	times::Array{Float64,1}
	discretization::Array{Float64,1}
	queue::Array{Float64,1}
	density::Array{Array{Float64,1},1}
	interval_lengths::Array{Float64,1}
	entrance::Junction
	exit::Junction
	ProductionLine(velocity, capacity, entrance, exit, queue, density::Function=x->0, cells=10^3) = begin
		flux(t,rho) = min(velocity * rho, capacity(t))
		len = exit - entrance
		Δx = len/cells
	
		discr = [ entrance + Δx/2 + n*Δx for n in 0:(cells-1) ]
		l = length(discr)
		borders = Array{Float64,1}(undef, length(discr)+1)
		borders[2:end-1] = [
			(discr[idx]+discr[idx+1])/2 
			for idx in 1:(length(discr)-1)
		]
		borders[1] = entrance
		borders[end] = exit
	
		interval_lengths = diff(borders)
	
		return new(
			flux, capacity, velocity, 
			[], discr, [queue], [density.(discr)], interval_lengths
		)
	end
end

# ╔═╡ 9c08ba10-3267-11eb-1309-3f25a24d65d0
struct ProductionNetwork <: Network
	junctions::Array{Junction,1}
	lines::Array{Line,1}
	max_timestep::Float64
	ProductionNetwork(junctions, lines) = begin
		t = maximum([maximum(line.interval_lengths)/line.velocity for line in lines])
		return new(junctions, lines, t)
	end
end

# ╔═╡ f62ca280-3268-11eb-1d2b-d7cb0a36a7bb
struct SourceJunction <: Junction
	exit::Line
	Flux
	SourceJunction(exit, Flux) = begin
		jun = new(exit, Flux)
		exit.entrance = jun
		return jun
	end
end

# ╔═╡ 4dc63a90-34af-11eb-2ffd-0d140c8a1504
function junction_flux(junction::SourceJunction, time_idx, time)
	return Dict(junction.exit => junction.Flux(time))
end

# ╔═╡ 93ca6c4e-34af-11eb-2d54-2364a39b0e72
struct InnerJunction <: Junction
	entrances::Array{Line,1}
	distribution::Dict{Line, Function}
	InnerJunction(entrances, distribution) = begin
		jun = new(entrances, distribution)
		for entrance in entrances
			entrance.exit = jun
		end
		for (exit, _) in distribution
			exit.entrance = jun
		end
		return jun
	end
end

# ╔═╡ 19643c62-34b0-11eb-074d-0356883e7ed2
function junction_flux(junction::InnerJunction, time_idx, time)
	total_influx = sum(
		[line.flux(time, line.density[time_idx]) for line in junction.entrances]
	)
	return Dict( line => A(time)*total_influx for (line, A) in junction.distribution)
end

# ╔═╡ 59779960-34b4-11eb-0664-953c5c5ecc9f
struct SinkJunction <: Junction
	entrance::Line
	SinkJunction(entrance) = begin
		jun = new(entrance)
		entrance.exit = jun
		return jun
	end
end

# ╔═╡ bbb94150-34b4-11eb-1c58-b9e0dec4f8c1
function junction_flux(junction::SinkJunction, time_idx, time)
	return Dict{Line, Float64}()
end

# ╔═╡ 5fb56360-34b5-11eb-0832-7ba19f79a05a
function post_influx(influx, capacity, queue, time_delta)
	if queue > 0
		return min(capacity, queue/time_delta + influx)
	else
		return min(influx, capacity)
	end
end

# ╔═╡ 53e34410-3269-11eb-13ed-0fdb12e71a14
function godunov_step!(
	line::ProductionLine, influx::Float64, time_idx,
)
	t0 = line.times[time_idx]
	time_delta = line.times[time_idx+1]-t0
	density = line.density[time_idx]
	queue = line.queue[time_idx]
	
	# calculate flux
	flux = Array{Float64,1}(undef, length(density)+1)
	flux[1] = post_influx(influx, line.capacity(t0), queue, time_delta)
	flux[2:end] = line.flux.(t0, density)
	
	increments = -diff(flux) # influx into a cell
	
	# update density and queue
	new_density = density .+ (time_delta ./ line.interval_lengths) .* increments
	line.density[time_idx+1] = new_density
	
	line.queue[time_idx+1] = queue + time_delta * (influx - flux[1])
end

# ╔═╡ b6ced280-3267-11eb-1979-b30e9214c173
function godunov_solve(network::ProductionNetwork, time_horizon)
		
	timesteps_required = ceil(time_horizon/network.max_timestep)
	timestep = time_horizon/timesteps_required
	number_of_timesteps = ceil(Int, time_horizon/timestep)
	times = [n*timestep for n in 0:number_of_timesteps]
	
	for line in network.lines
		append!(line.density, Array{Array{Float64,1},1}(undef, number_of_timesteps))
		append!(line.queue, Array{Float64,1}(undef, number_of_timesteps))
		line.times = times
	end
	
	for time_idx in 1:(number_of_timesteps)
		junc_flux = Dict{Junction, Dict{Line,Float64}}()
		for junction in network.junctions
			junc_flux[junction] = junction_flux(junction, time_idx, times[time_idx])
		end
		
		for line in network.lines
			influx = junc_flux[line.entrance][line]
			godunov_step!(line, influx, time_idx)
		end
	end
	
	return network.lines
end

# ╔═╡ b6559160-34b8-11eb-02e2-e76d36f46a79
md"# Exercise 3"

# ╔═╡ 61b6f860-34b8-11eb-255a-33acc897a55e
function density(line::ProductionLine, x, t)
	x_idx = max(searchsortedlast(line.discretization, x),1)
	t_idx = max(searchsortedlast(line.times, t),1)
	return line.density[t_idx][x_idx]
end

# ╔═╡ 0a557960-34b9-11eb-3ee4-0b8048f40a7c
function queue(line::ProductionLine, t)
	t_idx = max(searchsortedlast(line.times, t),1)
	return line.queue[t_idx]
end

# ╔═╡ 7ce74ede-34b9-11eb-30c2-b7e95de3e854
begin
	velocity = 1
	capacity(t) = 1
	queue0 = 0
	line0 = ProductionLine(velocity, capacity, 0, 1, queue0)
	source = SourceJunction(line0, t->1)
	sink = SinkJunction(line0)
	network = ProductionNetwork([source, sink], [line0])
	godunov_solve(network, 2)
end

# ╔═╡ 3dcc5770-34bc-11eb-3146-69718faf2fcc
plot(line0.discretization, line0.density[1001], ylim=(0,1))

# ╔═╡ Cell order:
# ╠═ccd42400-3266-11eb-3d2c-97c96caf3e78
# ╠═9863f500-3267-11eb-32e7-b5a82fb92dc9
# ╠═61bc4c92-3268-11eb-1197-7786485b47d1
# ╠═e31a36f0-3284-11eb-08a3-77c7fc3c47e6
# ╠═663ebce2-3267-11eb-1394-318abf53bd13
# ╠═9c08ba10-3267-11eb-1309-3f25a24d65d0
# ╠═f62ca280-3268-11eb-1d2b-d7cb0a36a7bb
# ╠═4dc63a90-34af-11eb-2ffd-0d140c8a1504
# ╠═93ca6c4e-34af-11eb-2d54-2364a39b0e72
# ╠═19643c62-34b0-11eb-074d-0356883e7ed2
# ╠═59779960-34b4-11eb-0664-953c5c5ecc9f
# ╠═bbb94150-34b4-11eb-1c58-b9e0dec4f8c1
# ╠═b6ced280-3267-11eb-1979-b30e9214c173
# ╠═53e34410-3269-11eb-13ed-0fdb12e71a14
# ╠═5fb56360-34b5-11eb-0832-7ba19f79a05a
# ╟─b6559160-34b8-11eb-02e2-e76d36f46a79
# ╠═61b6f860-34b8-11eb-255a-33acc897a55e
# ╠═0a557960-34b9-11eb-3ee4-0b8048f40a7c
# ╠═7ce74ede-34b9-11eb-30c2-b7e95de3e854
# ╠═3dcc5770-34bc-11eb-3146-69718faf2fcc
