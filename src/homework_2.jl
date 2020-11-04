### A Pluto.jl notebook ###
# v0.12.6

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d456a4d0-1c33-11eb-3c70-599cd64c2c83
begin
	using Printf: @sprintf
	using Zygote
	using Plots: plot, plot!
	using Roots: find_zero
	using PlutoUI
end

# ╔═╡ 1f4955c0-1b9c-11eb-2e53-c10ebcb1ff07
md"# Exercise 1"

# ╔═╡ 216e96c0-1c4c-11eb-0ae0-d7af7fc2f12e
function rho_star_from_flux(flux, flux_extremum)
	dflux = flux'
	return (rho_left, rho_right) -> begin
		s_negative = ((flux(rho_left)-flux(rho_right))<0) ⊻ ((rho_left-rho_right)<0)
		
		f_left = dflux(rho_left)
		f_right = dflux(rho_right)
		if f_left > 0 
			if f_right > 0
				return rho_left
			elseif s_negative
				return rho_right
			else
				return rho_left
			end
		else
			if f_right > 0
				return flux_extremum
			else
				return rho_right
			end
		end
	end
end

# ╔═╡ 26307230-1dd9-11eb-2d5c-691039ab0dde
abstract type Junction end

# ╔═╡ 66da2a4e-1db2-11eb-2c28-a3b282c77228
mutable struct Road
	flux
	flux_max::Float64
	rho_star
	discretization::Array{Float64,1}
	density::Array{Array{Float64,1},1}
	interval_lengths::Array{Float64,1}
	entrance::Junction
	exit::Junction
	Road(flux, entrance, exit, density::Function, cells=10^3) = begin
		flux_max = find_zero(flux', 0)
		rho_star = rho_star_from_flux(flux, flux_max)
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
			flux, flux_max, rho_star, discr, [density.(discr)], interval_lengths
		)
	end
end

# ╔═╡ 692139b0-1dea-11eb-0696-278fb897307d
struct Network
	junctions::Array{Junction,1}
	roads::Array{Road,1}
	max_timestep::Float64
	Network(junctions, roads) = begin
		running_min = Inf
		for road in roads
			u_den = unique(road.density[1])
			abs_derivatives = abs.(road.flux'.(u_den))
			norm_dflux = 2*maximum(abs_derivatives)
			running_min = min(running_min, minimum(road.interval_lengths)/norm_dflux)
		end
		return new(junctions, roads, running_min)
	end
end

# ╔═╡ b7bff5b0-1de1-11eb-38f4-8be3fcb3424d
struct JunctionFlux
	entrances::Dict{Road, Float64}
	exits::Dict{Road,Float64}
end

# ╔═╡ 2ef05130-1df0-11eb-38bf-9ba9c9fa4590
struct SourceJunction <: Junction
	exit::Road
	constantFlux
end

# ╔═╡ 9732d6a0-1df0-11eb-17ef-b974a53594d3
struct SinkJunction <: Junction
	entrance::Road
	constantFlux
end

# ╔═╡ 04cdeba0-1df1-11eb-018d-293f4a6bcf1a
function junction_flux(junction::SinkJunction, time_idx)
	entrance_flux = Dict(junction.entrance => junction.constantFlux)
	return JunctionFlux(entrance_flux, Dict{Road,Float64}())
end

# ╔═╡ 9f89bca0-1df1-11eb-0545-f759d4ab800e
function junction_flux(junction::SourceJunction, time_idx)
	exit_flux = Dict(junction.exit => junction.constantFlux)
	return JunctionFlux(Dict{Road,Float64}(), exit_flux)
end

# ╔═╡ 52974242-1c53-11eb-1002-f125de4500c0
md"## Beispiel"

# ╔═╡ e6a74fc0-1c53-11eb-3bd2-d3eab118a1bd
md"### $\rho_L$"

# ╔═╡ 5fb2f8c0-1c53-11eb-35e0-2b61511a6dbf
@bind ρ_L Slider(0.05:0.05:1, default = 0.2, show_value=true)

# ╔═╡ db30bdc0-1c53-11eb-028d-e108efac2e3d
md"### $\rho_R$"

# ╔═╡ ca4820c0-1c53-11eb-39ac-170b702e43b1
@bind ρ_R Slider(0.05:0.05:1, default = 0.5, show_value=true)

# ╔═╡ 83e3ff40-1c54-11eb-2810-a5f306cf9557
md"### Time Horizon"

# ╔═╡ 8c81f5d0-1c54-11eb-00f3-2fbd04ded608
@bind time_horizon Slider(0.2:0.2:6, default=1, show_value=true)

# ╔═╡ 3abd2ea0-1df2-11eb-110d-95d6032829e9
begin
	flux(ρ) = ρ*(1-ρ)
	
	density_0(x) = x<0 ? ρ_L : ρ_R 
	
	road = Road(flux, -5, 5, density_0)
	source_junction = SourceJunction(road, flux(ρ_L))
	sink_junction = SinkJunction(road, flux(ρ_R))
	road.entrance = source_junction
	road.exit = sink_junction
	network = Network([source_junction, sink_junction], [road])
end

# ╔═╡ b6d1fd20-1dea-11eb-2e8d-71fe7fb602b7
function godunov_step(
	road::Road, influx::Float64, outflux::Float64, density, time_delta
)
	flux = Array{Float64,1}(undef, length(density)+1) # n cells have n+1 borders
	flux[1] = influx
	flux[end] = outflux
	flux[2:end-1] = [
		road.flux(road.rho_star(x,y)) 
		for (x,y) in zip(density[1:end-1], density[2:end])
	]
	
	increments = -diff(flux) # influx into a cell
	return density .+ (time_delta ./ road.interval_lengths) .* increments
end

# ╔═╡ b921f3b0-1e02-11eb-2dae-b5693f16837f
begin
	discr = road.discretization
	density = road.density
end

# ╔═╡ 9be3f732-1c72-11eb-0e6f-a1d96e46797a
md"# Exercise2"

# ╔═╡ 3f95ccb0-1dda-11eb-245d-03b638c1c4eb
function demand(road::Road, time_idx)
	rho = road.density[time_idx][end]
	if rho < road.flux_max
		return road.flux(rho)
	else
		return road.flux(road.flux_max)
	end
end

# ╔═╡ 6abbd4b0-1ddb-11eb-027a-7762903061f3
function supply(road::Road, time_idx)
	rho = road.density[time_idx][1]
	if rho < road.flux_max
		return road.flux(road.flux_max)
	else
		return road.flux(rho)
	end
end

# ╔═╡ 7cc675c0-1dd6-11eb-240e-3193d1ace5ed
struct OneTwoJunction <: Junction
	entrances::Array{Road,1}
	exits::Array{Road,1}
	distribution_matrix::Array{Float64,2}
end

# ╔═╡ e457cc40-1dd9-11eb-1388-75ea1e06bca7
function junction_flux(junction::OneTwoJunction, time_idx)
	d = [demand(road, time_idx) for road in junction.entrances]
	s = [supply(road, time_idx) for road in junction.exits]
	
	entrance_flux = Array{Float64,1}(undef, length(junction.entrances))
	
	#abuse one-to-two special case knowledge for solution
	entrance_flux[1] = min(d[1], minimum(s ./ junction.distribution_matrix[:,1]))
	exit_flux = junction.distribution_matrix[:,1] .* entrance_flux[1]
	##
	
	return JunctionFlux(entrance_flux, exit_flux)
end

# ╔═╡ efbd74b0-1de1-11eb-2875-d5e772b9625d
struct TwoOneJunction <: Junction
	entrances::Array{Road,1}
	exits::Array{Road,1}
	right_of_way::Float64
end

# ╔═╡ 431e84a0-1de2-11eb-33ba-b9abe149f86f
function junction_flux(junction::TwoOneJunction, time_idx)
	d = [demand(road, time_idx) for road in junction.entrances]
	s = [supply(road, time_idx) for road in junction.exits]
	
	# abuse two-to-one special case knowledge for solution (Right of way priority)
	q = junction.right_of_way
	maximal_flow = min(d[1]+d[2], s[1])
	
	entrance_flux[1] = min(d[1], q/(1-q) * d[2], q*maximal_flow)
	entrance_flux[2] = (1-q)/q * entrance_flux[1]
	
	exit_flux = sum(entrance_flux)
	##
	
	return JunctionFlux(entrance_flux, exit_flux)
end

# ╔═╡ 047bcba0-1dff-11eb-145e-ddc4933e8826
function godunov_solve(network::Network, time_horizon)
		
	timesteps_required = ceil(time_horizon/network.max_timestep)
	timestep = time_horizon/timesteps_required
	number_of_timesteps = ceil(Int, time_horizon/timestep)
	
	for road in network.roads
		append!(road.density, Array{Array{Float64,1},1}(undef, number_of_timesteps))
	end
	
	for time_idx in 1:(number_of_timesteps)
		junc_flux = Dict{Junction, JunctionFlux}()
		for junction in network.junctions
			junc_flux[junction] = junction_flux(junction, time_idx)
		end
		
		for road in network.roads
			influx = junc_flux[road.entrance].exits[road]
			outflux = junc_flux[road.exit].entrances[road]
			road.density[time_idx+1]=godunov_step(
				road, influx, outflux, road.density[time_idx], timestep
			)
		end
	end
	
	times = [n*timestep for n in 0:number_of_timesteps]
	return times, network.roads
end

# ╔═╡ 46091370-1e7c-11eb-0943-ad3fb2982453
(times, roads) = godunov_solve(network, time_horizon)

# ╔═╡ 6aaa0640-1c6e-11eb-31dd-d7223b7df72a
idx_of_t(t) = findfirst(x->x>=t,times)

# ╔═╡ 11073d10-1c5f-11eb-078b-714e16d80478
begin
	p = plot()
	for t in 0:time_horizon/3-0.00001:time_horizon
		plot!(
			discr, density[idx_of_t(t)];
			linetype=:steppost, label="t=$(@sprintf "%0.2f" t)" 
		)
	end
	p
end

# ╔═╡ 7c5be130-1de3-11eb-14c4-21e5ef54a48a
md"# Exercise 3"

# ╔═╡ Cell order:
# ╠═d456a4d0-1c33-11eb-3c70-599cd64c2c83
# ╟─1f4955c0-1b9c-11eb-2e53-c10ebcb1ff07
# ╠═216e96c0-1c4c-11eb-0ae0-d7af7fc2f12e
# ╠═26307230-1dd9-11eb-2d5c-691039ab0dde
# ╠═66da2a4e-1db2-11eb-2c28-a3b282c77228
# ╠═692139b0-1dea-11eb-0696-278fb897307d
# ╠═b7bff5b0-1de1-11eb-38f4-8be3fcb3424d
# ╠═2ef05130-1df0-11eb-38bf-9ba9c9fa4590
# ╠═9732d6a0-1df0-11eb-17ef-b974a53594d3
# ╠═04cdeba0-1df1-11eb-018d-293f4a6bcf1a
# ╠═9f89bca0-1df1-11eb-0545-f759d4ab800e
# ╠═047bcba0-1dff-11eb-145e-ddc4933e8826
# ╠═b6d1fd20-1dea-11eb-2e8d-71fe7fb602b7
# ╟─52974242-1c53-11eb-1002-f125de4500c0
# ╟─e6a74fc0-1c53-11eb-3bd2-d3eab118a1bd
# ╟─5fb2f8c0-1c53-11eb-35e0-2b61511a6dbf
# ╟─db30bdc0-1c53-11eb-028d-e108efac2e3d
# ╠═ca4820c0-1c53-11eb-39ac-170b702e43b1
# ╟─83e3ff40-1c54-11eb-2810-a5f306cf9557
# ╟─8c81f5d0-1c54-11eb-00f3-2fbd04ded608
# ╠═3abd2ea0-1df2-11eb-110d-95d6032829e9
# ╠═46091370-1e7c-11eb-0943-ad3fb2982453
# ╠═b921f3b0-1e02-11eb-2dae-b5693f16837f
# ╠═6aaa0640-1c6e-11eb-31dd-d7223b7df72a
# ╠═11073d10-1c5f-11eb-078b-714e16d80478
# ╟─9be3f732-1c72-11eb-0e6f-a1d96e46797a
# ╠═3f95ccb0-1dda-11eb-245d-03b638c1c4eb
# ╠═6abbd4b0-1ddb-11eb-027a-7762903061f3
# ╠═7cc675c0-1dd6-11eb-240e-3193d1ace5ed
# ╠═e457cc40-1dd9-11eb-1388-75ea1e06bca7
# ╠═efbd74b0-1de1-11eb-2875-d5e772b9625d
# ╠═431e84a0-1de2-11eb-33ba-b9abe149f86f
# ╟─7c5be130-1de3-11eb-14c4-21e5ef54a48a
