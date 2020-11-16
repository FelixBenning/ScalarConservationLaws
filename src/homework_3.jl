### A Pluto.jl notebook ###
# v0.12.10

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

# ╔═╡ f486b5f0-283c-11eb-2290-83e6460d3546
begin
	using PlutoUI:Slider
end

# ╔═╡ 62d0ee10-2679-11eb-25f8-3fda2717e2fd
md"# Exercise 1"

# ╔═╡ 1bbc7c70-2687-11eb-3c08-6504df16ad41
function arrival_times(f_A, eps, N)
	times = Array{Float64,1}(undef, N+1)
	times[1] = 0
	for idx in 1:N
		times[idx+1] = times[idx] + eps/f_A(times[idx])  
	end
	return times
end

# ╔═╡ 902052ce-2687-11eb-28cc-5d5c30074f59
function arrival_grid(tau_A::Union{Array,UnitRange}, eps, mu, N, M, K)
	grid = Array{Float64, 2}(undef, K*M+1, N+1)
	grid[1,:] = tau_A

	for processor_idx in 1:K*M
		m = floor(Int, (processor_idx-1)/K)
		grid[processor_idx+1,1] = grid[processor_idx,1] + 1/mu[m+1]
		for item_idx in 2:N+1
			grid[processor_idx+1, item_idx] = max(
				grid[processor_idx, item_idx], 
				grid[processor_idx+1, item_idx-1]
			) + eps/mu[m+1]
		end
	end
	return grid
end

# ╔═╡ d8a81030-274d-11eb-27f0-7107380751f7
function arrival_grid(f_A::Function, eps, mu, N, M, K)
	return arrival_grid(arrival_times(f_A, eps,N), eps, mu, N, M, K)
end

# ╔═╡ d8009be2-2769-11eb-28f2-a71305a3ebd9
md"### Test with 6.3"

# ╔═╡ fe586f72-2769-11eb-2c61-fb71e75cb60e
begin
	N_1 = 3
	M_1 = 2
	τ_A = 0:N_1
	μ = [1,0.5,1]
	arrival_grid(τ_A, 1, μ, N_1, M_1, 1)
end

# ╔═╡ 779611f0-283a-11eb-1e09-b9cc51a6f71e
md"# Exercise 2"

# ╔═╡ e1703c6e-283c-11eb-1dc3-b7ce3491cfe3
md"#### N"

# ╔═╡ edead0f0-283c-11eb-0403-0b46d7c8dc4c
@bind N Slider(1:100, default = 10, show_value=true)

# ╔═╡ 320ea9f0-283d-11eb-37e5-d5eca83ae37c
md"#### K"

# ╔═╡ 36496d6e-283d-11eb-0352-0f2b2518dbe5
@bind K Slider(1:100, default = 10, show_value=true)

# ╔═╡ 2a0dec82-2837-11eb-3dc6-3d3724864e7c
begin
	v = [2,1,1/2]
	f_A(t) = 3/2
	T = 5
	eps = 1/K
	ρ = density_from_arrivals(arrival_grid(f_A, 1/K, v, N, 3, K))
end

# ╔═╡ 070273f2-276a-11eb-227d-6b3259d00b97
function density_from_arrivals(arrival_grid)
	return (x, t) -> begin #ρ_eps density function
		k = floor(Int, x/eps) # processor index
		u_k = searchsortedfirst(time_grid[k,:], t) # index is number of goods passed
		u_k_plus_1 = searchsortedfirst(time_grid[k+1,:], t)
		return u_k - u_k_plus_1 # work in progress at (virtual) processor k
	end
end

# ╔═╡ Cell order:
# ╟─62d0ee10-2679-11eb-25f8-3fda2717e2fd
# ╠═f486b5f0-283c-11eb-2290-83e6460d3546
# ╠═1bbc7c70-2687-11eb-3c08-6504df16ad41
# ╠═902052ce-2687-11eb-28cc-5d5c30074f59
# ╠═d8a81030-274d-11eb-27f0-7107380751f7
# ╟─d8009be2-2769-11eb-28f2-a71305a3ebd9
# ╠═fe586f72-2769-11eb-2c61-fb71e75cb60e
# ╟─779611f0-283a-11eb-1e09-b9cc51a6f71e
# ╠═070273f2-276a-11eb-227d-6b3259d00b97
# ╟─e1703c6e-283c-11eb-1dc3-b7ce3491cfe3
# ╟─edead0f0-283c-11eb-0403-0b46d7c8dc4c
# ╟─320ea9f0-283d-11eb-37e5-d5eca83ae37c
# ╟─36496d6e-283d-11eb-0352-0f2b2518dbe5
# ╠═2a0dec82-2837-11eb-3dc6-3d3724864e7c
