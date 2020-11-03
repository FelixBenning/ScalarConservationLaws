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

# ╔═╡ 41175b30-1c3b-11eb-10d5-91414cfbb583
function _godunov_solve(
	flux, start_flux, stop_flux, discr, density_0, interval_lengths, times, rho_star
)
	density = Array{Float64, 2}(undef, length(density_0), length(times))
	density[:,1] = density_0
	for (idx, delta_t) in enumerate(diff(times))
		den = density[:,idx]
		
		# flux on the interval borders
		b_fl = Array{Float64,1}(undef, length(density_0)+1)
		b_fl[2:end-1]=[flux(rho_star(x,y)) for (x,y) in zip(den[1:end-1], den[2:end])]
		b_fl[1] = start_flux(times[idx])
		b_fl[end] = stop_flux(times[idx])
		
		incr = -diff(b_fl) # influx into a cell
		density[:,idx+1] = den .+ delta_t./interval_lengths .* incr
	end
	return (times, discr, density)
end

# ╔═╡ 216e96c0-1c4c-11eb-0ae0-d7af7fc2f12e
function rho_star_from_flux(flux, roothint=0)
	dflux = flux'
	dflux_inv_0 = find_zero(dflux, roothint)

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
				return dflux_inv_0
			else
				return rho_right
			end
		end
	end
end

# ╔═╡ 430325e0-1b9c-11eb-26f7-93813a1b9698
"""
flux: a (concave or convex, twice differentiable) flux function which accepts a density as an argument
start_flux: flux(ρ(start,t)) = start_flux(t)
stop_flux: flux(ρ(stop, t)) = stop_flux(t)
start, stop: ends of the interval
discr: discretization points matching density_0
density_0: the density at time 0 evaluated at the discretization points
time_horizon: calculate flow until time_horizon
"""
function discrete_godunov_solve(
	flux, start_flux, stop_flux, start, stop, discr, density_0, time_horizon
)
	borders = Array{Float64,1}(undef, length(discr)+1)
	borders[2:end-1] = [(discr[idx]+discr[idx+1])/2 for idx in 1:(length(discr)-1)]
	borders[1] = start
	borders[end] = stop
	interval_lengths = diff(borders)
	
	u_den = unique(density_0)
	abs_derivatives = abs.(flux'.(u_den))
	max_time_step = minimum(interval_lengths)/(2*maximum(abs_derivatives))
	rho_star = rho_star_from_flux(flux, u_den[argmin(abs_derivatives)]) 
	
	timesteps_required = ceil(time_horizon/max_time_step)
	time_step = time_horizon/timesteps_required
	times = [n*time_step for n in 0:ceil(time_horizon/time_step)]
	return _godunov_solve(
		flux, start_flux, stop_flux, discr, density_0, 
		interval_lengths, times, rho_star
	)
end

# ╔═╡ 60952d50-1c51-11eb-09ed-115f0cb1f8d7
"""
flux: a (concave or convex, twice differentiable) flux function which accepts a density as an argument
start_flux: flux(ρ(start,t)) = start_flux(t)
stop_flux: flux(ρ(stop, t)) = stop_flux(t)
start, stop: ends of the interval
discr: discretization points matching density_0
density_0: the density function at time 0
time_horizon: calculate flow until time_horizon
cells: number of cells to approximate density (default 10 000)
"""
function godunov_solve(
	flux, start_flux, stop_flux, start, stop, density_0, time_horizon, cells = 10^3
)
	L = stop - start
	Δx = L / cells
	discr = [ start + Δx/2 + n*Δx for n in 0:(cells-1)]
	
	return discrete_godunov_solve(
		flux, start_flux, stop_flux, start, stop, 
		discr, density_0.(discr), time_horizon
	)
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

# ╔═╡ 2441afc0-1c53-11eb-2f42-215a417f4dde
begin
	flux(ρ) = ρ*(1-ρ)
	start_flux(t) = flux(ρ_L)
	stop_flux(t) = flux(ρ_R)
	start = -5
	stop = 5
	
	density_0(x) = x<0 ? ρ_L : ρ_R 
	
	times, discr, density = godunov_solve(
		flux, start_flux, stop_flux, start, stop, density_0, time_horizon)
end

# ╔═╡ 6aaa0640-1c6e-11eb-31dd-d7223b7df72a
idx_of_t(t) = findfirst(x->x>=t,times)

# ╔═╡ 11073d10-1c5f-11eb-078b-714e16d80478
begin
	p = plot()
	for t in 0:time_horizon/3-0.00001:time_horizon
		plot!(
			discr, density[:,idx_of_t(t)];
			linetype=:steppost, label="t=$(@sprintf "%0.2f" t)" 
		)
	end
	p
end

# ╔═╡ 9be3f732-1c72-11eb-0e6f-a1d96e46797a
md"# Exercise2"

# ╔═╡ 26307230-1dd9-11eb-2d5c-691039ab0dde
abstract type Junction end

# ╔═╡ 66da2a4e-1db2-11eb-2c28-a3b282c77228
struct Road
	flux
	flux_max::Float64
	discretization::Array{Float64,1}
	interval_lengths::Array{Float64,1}
	density::Array{Float64,2}
	entrance::Junction
	exit::Junction
end

# ╔═╡ 3f95ccb0-1dda-11eb-245d-03b638c1c4eb
function demand(road::Road, time_idx)
	rho = road.density[end,time_idx]
	if rho < road.flux_max
		return road.flux(rho)
	else
		return road.flux(road.flux_max)
	end
end

# ╔═╡ 6abbd4b0-1ddb-11eb-027a-7762903061f3
function supply(road::Road, time_idx)
	rho = road.density[1,time_idx]
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

# ╔═╡ b7bff5b0-1de1-11eb-38f4-8be3fcb3424d
struct JunctionFlux
	entrances::Array{Float64,1}
	exits::Array{Float64,1}
end

# ╔═╡ e457cc40-1dd9-11eb-1388-75ea1e06bca7
function junction_flux(junction::OneTwoJunction, time)
	d = [demand(road, time) for road in junction.entrances]
	s = [supply(road, time) for road in junction.exits]
	
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
function junction_flux(junction::TwoOneJunction, time)
	d = [demand(road, time) for road in junction.entrances]
	s = [supply(road, time) for road in junction.exits]
	
	# abuse two-to-one special case knowledge for solution (Right of way priority)
	q = junction.right_of_way
	maximal_flow = min(d[1]+d[2], s[1])
	
	entrance_flux[1] = min(d[1], q/(1-q) * d[2], q*maximal_flow)
	entrance_flux[2] = (1-q)/q * entrance_flux[1]
	
	exit_flux = sum(entrance_flux)
	##
	
	return JunctionFlux(entrance_flux, exit_flux)
end

# ╔═╡ 7c5be130-1de3-11eb-14c4-21e5ef54a48a


# ╔═╡ Cell order:
# ╠═d456a4d0-1c33-11eb-3c70-599cd64c2c83
# ╟─1f4955c0-1b9c-11eb-2e53-c10ebcb1ff07
# ╠═41175b30-1c3b-11eb-10d5-91414cfbb583
# ╠═216e96c0-1c4c-11eb-0ae0-d7af7fc2f12e
# ╠═430325e0-1b9c-11eb-26f7-93813a1b9698
# ╠═60952d50-1c51-11eb-09ed-115f0cb1f8d7
# ╟─52974242-1c53-11eb-1002-f125de4500c0
# ╟─e6a74fc0-1c53-11eb-3bd2-d3eab118a1bd
# ╟─5fb2f8c0-1c53-11eb-35e0-2b61511a6dbf
# ╟─db30bdc0-1c53-11eb-028d-e108efac2e3d
# ╠═ca4820c0-1c53-11eb-39ac-170b702e43b1
# ╟─83e3ff40-1c54-11eb-2810-a5f306cf9557
# ╟─8c81f5d0-1c54-11eb-00f3-2fbd04ded608
# ╠═2441afc0-1c53-11eb-2f42-215a417f4dde
# ╠═6aaa0640-1c6e-11eb-31dd-d7223b7df72a
# ╠═11073d10-1c5f-11eb-078b-714e16d80478
# ╟─9be3f732-1c72-11eb-0e6f-a1d96e46797a
# ╠═26307230-1dd9-11eb-2d5c-691039ab0dde
# ╠═3f95ccb0-1dda-11eb-245d-03b638c1c4eb
# ╠═6abbd4b0-1ddb-11eb-027a-7762903061f3
# ╠═66da2a4e-1db2-11eb-2c28-a3b282c77228
# ╠═7cc675c0-1dd6-11eb-240e-3193d1ace5ed
# ╠═b7bff5b0-1de1-11eb-38f4-8be3fcb3424d
# ╠═e457cc40-1dd9-11eb-1388-75ea1e06bca7
# ╠═efbd74b0-1de1-11eb-2875-d5e772b9625d
# ╠═431e84a0-1de2-11eb-33ba-b9abe149f86f
# ╠═7c5be130-1de3-11eb-14c4-21e5ef54a48a
