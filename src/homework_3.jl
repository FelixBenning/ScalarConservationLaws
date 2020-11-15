### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 62d0ee10-2679-11eb-25f8-3fda2717e2fd
md"# Exercise 1"

# ╔═╡ 1bbc7c70-2687-11eb-3c08-6504df16ad41
function arrival_times(f_A, eps, N)
	times = Array{Float64,1}(undef, N)
	times[1] = 0
	for idx in 1:N-1
		times[idx+1] = times[idx] + eps/f_A(times[idx])  
	end
	return times
end

# ╔═╡ 902052ce-2687-11eb-28cc-5d5c30074f59
function arrival_grid(tau_A::Union{Array,UnitRange}, eps, mu, N, M, K)
	grid = Array{Float64, 2}(undef, K*M+1, N+1)
	grid[1,:] = tau_A

	for processor_idx in 2:K*M+1
		m = floor(Int, processor_idx/K)
		grid[processor_idx,1] = grid[processor_idx-1,1] + 1/mu[m]
		for item_idx in 2:N+1
			grid[processor_idx, item_idx] = max(
				grid[processor_idx-1, item_idx], 
				grid[processor_idx, item_idx-1]
			) + eps/mu[m]
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
	N = 3
	M = 2
	τ_A = 0:N
	μ = [1,0.5,1]
	arrival_grid(τ_A, 1, μ, N, M, 1)
end

# ╔═╡ 070273f2-276a-11eb-227d-6b3259d00b97


# ╔═╡ Cell order:
# ╟─62d0ee10-2679-11eb-25f8-3fda2717e2fd
# ╠═1bbc7c70-2687-11eb-3c08-6504df16ad41
# ╠═902052ce-2687-11eb-28cc-5d5c30074f59
# ╠═d8a81030-274d-11eb-27f0-7107380751f7
# ╟─d8009be2-2769-11eb-28f2-a71305a3ebd9
# ╠═fe586f72-2769-11eb-2c61-fb71e75cb60e
# ╠═070273f2-276a-11eb-227d-6b3259d00b97
