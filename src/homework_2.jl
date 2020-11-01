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
	using Zygote
	using Plots: plot
	using Roots: find_zero
	using PlutoUI
end

# ╔═╡ 1f4955c0-1b9c-11eb-2e53-c10ebcb1ff07
md"# Exercise 1"

# ╔═╡ 41175b30-1c3b-11eb-10d5-91414cfbb583
function _godunov_solve(
	flux, start_flux, stop_flux, discr, density_0, interval_lengths, times, rho_star
)
	time_diffs = diff(times)
	

	density = Array{Float64, 2}(undef, length(density_0), length(times))
	density[:,1] = density_0
	for (idx, delta_t) in enumerate(time_diffs)
		den = density[:,idx]
		b_fl = Array{Float64,1}(undef, length(density_0)+1)
		b_fl[2:end-1]=[flux(rho_star(x,y)) for (x,y) in zip(den[1:end-1], den[2:end])]
		
		b_fl[1] = start_flux(times[idx])
		b_fl[end] = stop_flux(times[idx])
		
		incr = -diff(b_fl)
		density[:,idx+1] = den .+ delta_t./interval_lengths .* incr
	end
	return (times, discr, density)
end

# ╔═╡ bce334d0-1c6a-11eb-2ea5-454581da7761
Array{Float64,2}(undef, 2,3)

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
	interval_borders = [(discr[idx]+discr[idx+1])/2 for idx in 1:(length(discr)-1)]
	prepend!(interval_borders, start)
	append!(interval_borders, stop)
	interval_lengths = diff(interval_borders)
	
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
@bind time_horizon Slider(0.2:0.2:3, default=1, show_value=true)

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

# ╔═╡ 11073d10-1c5f-11eb-078b-714e16d80478
plot(discr, density[:,end], linetype=:steppost)

# ╔═╡ b43e09c2-1c53-11eb-1d5a-ef39a0013ff9
false ? 5 : 2

# ╔═╡ 72c664d0-1c2e-11eb-1ddc-b1282c7e8e32
a =[1.0,2,3]

# ╔═╡ 20161530-1c58-11eb-222a-2d4473330f47
typeof(a)

# ╔═╡ df8ef550-1c33-11eb-3c99-7fde08fab8b1
Δx = 2

# ╔═╡ e35cb870-1c2e-11eb-314c-b1395544a443
push!([1,2], [1,2][end] + 0)

# ╔═╡ Cell order:
# ╟─1f4955c0-1b9c-11eb-2e53-c10ebcb1ff07
# ╠═d456a4d0-1c33-11eb-3c70-599cd64c2c83
# ╠═41175b30-1c3b-11eb-10d5-91414cfbb583
# ╠═bce334d0-1c6a-11eb-2ea5-454581da7761
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
# ╠═11073d10-1c5f-11eb-078b-714e16d80478
# ╠═b43e09c2-1c53-11eb-1d5a-ef39a0013ff9
# ╠═72c664d0-1c2e-11eb-1ddc-b1282c7e8e32
# ╠═20161530-1c58-11eb-222a-2d4473330f47
# ╠═df8ef550-1c33-11eb-3c99-7fde08fab8b1
# ╠═e35cb870-1c2e-11eb-314c-b1395544a443
