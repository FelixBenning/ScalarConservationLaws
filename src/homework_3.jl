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
	using Plots:plot,plot!
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
function arrival_grid(tau_A::Union{Array,UnitRange}, eps, mu, T, N, M, K)
	grid = Array{Float64, 2}(undef, K*M+1, N+1)
	grid[1,:] = tau_A

	for processor_idx in 1:K*M
		m = floor(Int, (processor_idx-1)/K)
		grid[processor_idx+1,1] = grid[processor_idx,1] + 1/mu[m+1]
		for item_idx in 2:N+1
			grid[processor_idx+1, item_idx] = max(
				grid[processor_idx, item_idx] + eps*T[m+1],
				grid[processor_idx+1, item_idx-1] + eps/mu[m+1]
			)
		end
	end
	return grid
end

# ╔═╡ d8a81030-274d-11eb-27f0-7107380751f7
function arrival_grid(f_A::Function, eps, mu, T, N, M, K)
	return arrival_grid(arrival_times(f_A, eps,N), eps, mu, T, N, M, K)
end

# ╔═╡ d8009be2-2769-11eb-28f2-a71305a3ebd9
md"### Test with 6.3"

# ╔═╡ fe586f72-2769-11eb-2c61-fb71e75cb60e
begin
	N_1 = 3
	M_1 = 2
	τ_A = 0:N_1
	μ = [1,0.5,1]
	arrival_grid(τ_A, 1, μ, 1 ./μ, N_1, M_1, 1)
end

# ╔═╡ 779611f0-283a-11eb-1e09-b9cc51a6f71e
md"# Exercise 2"

# ╔═╡ 070273f2-276a-11eb-227d-6b3259d00b97
function density_from_arrivals(arrival_grid, eps)
	return (x, t) -> begin #ρ_eps density function
		k = floor(Int, x/eps)+1 # processor index
		u_k = searchsortedfirst(arrival_grid[k,:], t) #index is number of goods passed
		u_k_plus_1 = searchsortedfirst(arrival_grid[k+1,:], t)
		return u_k - u_k_plus_1 # work in progress at (virtual) processor k
	end
end

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
	M = 3
	mu = [2, 1, 2]
	v = [2,1,1/2]
	f_A(t) = 3/2
	eps = 1/K
	ρ = density_from_arrivals(arrival_grid(f_A, eps, mu, 1 ./v, N, M, K), eps)
	
	T = 5
	Δt = 2*eps/30
	times = 0:Δt:T
	queues = [ρ.(x,times) for x in 0:M-1]
end

# ╔═╡ c79f9e30-2855-11eb-2e82-2b9f4b94caf8
begin
	plots = [plot(times, q, label="q$(n)") for (n,q) in enumerate(queues)]
	plot(plots[1],plots[2],plots[3], layout=(1,3))
end

# ╔═╡ 9ccba9a0-285b-11eb-36da-b9b2eaad4a90
md"# Exercise 3"

# ╔═╡ a818baa0-285b-11eb-1766-99959ed836c0
begin
	queue_fun = [t->0, t->t/2*(t>1/2), t->0]
	density_fun = (x,t) -> begin
		result = 3/4*(0<x<2*t)*(0<x<1) # density on (0,1)
		result += (1<x<(t+1/2))*(1<x<2) # density on (1,2)
		result += 2*(2<x<(t/2+5/4))*(2<x<3) # density on (2,3)
		return result
	end
end

# ╔═╡ 404d9450-2c1a-11eb-1837-4512f56f93ed
md"#### (a)"

# ╔═╡ e6413a10-2bed-11eb-00f6-5d88ebe4a8be
function L1_error(T, delta_t, K, queue_1, density_1, queue_2, density_2)
	times = 0:delta_t:T

	M = length(queue_1)
	Δx = 1/K
	inner_locations = filter(x->!isapprox(x%1,0), 0:Δx:M) # filter out integers
	
	error = 0
	for t in times
		q_err = sum([abs(q1(t)-q2(t)) for (q1,q2) in zip(queue_1, queue_2)])
		d_err = sum([abs(density_1(x,t)-density_2(x,t)) for x in inner_locations])
		error += (q_err+d_err*Δx)*delta_t
	end
	return error
end

# ╔═╡ 26fc5480-2c0d-11eb-1aa8-f97eb83da977
filter(x->!isapprox(x%1,0.0), 0:1/10:3)

# ╔═╡ 98d3d782-2c09-11eb-0875-79be1af36267
@bind N_2 Slider(10:10:100, default=20, show_value=true)

# ╔═╡ 2c206df0-2c0a-11eb-1a9f-a10b0b8ee4c7
begin
	emp_dens(k) = density_from_arrivals(
		arrival_grid(f_A, 1/k, mu, 1 ./v, N_2, M, k), 
		1/k
	)
	Ks = 10:10:100
	del_t(k) = 2/(k*30)
	error = [
		L1_error(
			T, del_t(K1), K1, 
			[t->emp_dens(K1)(x,t) for x in 0:M-1], emp_dens(K1),
			queue_fun, density_fun
		)
		for K1 in Ks
	]
	plot(Ks, error)
end

# ╔═╡ 05f9e390-2c0a-11eb-1ce5-f1d4489d9273
md"#### (b)"

# ╔═╡ 956e8b70-2c14-11eb-0268-3b5827f67b7f
function integral(T, delta_t, K, queue, density, test_fun)
	times = 0:delta_t:T

	M = length(queue)
	Δx = 1/K
	inner_locations = filter(x->!isapprox(x%1,0), 0:Δx:M) # filter out integers
	
	error = 0
	for t in times
		q_err = sum([q(t)*test_fun(m-1,t) for (m,q) in enumerate(queue)])
		d_err = sum([density(x,t)*test_fun(x,t) for x in inner_locations])
		error += (q_err+d_err*Δx)*delta_t
	end
	return error
end

# ╔═╡ b7703960-2c16-11eb-3e0f-f3f194aa6f6a
@bind alpha_x Slider(0.4:0.1:4, default=2, show_value=true)

# ╔═╡ d4a9e9e0-2c16-11eb-1f56-8159f918e644
@bind alpha_t Slider(0.4:0.1:4, default=2, show_value=true)

# ╔═╡ 8d024be0-2c17-11eb-1e1f-41b299015af9
@bind x_c Slider(0.1:0.1:(3-1/alpha_x), default = 0.1, show_value=true)

# ╔═╡ f73e3050-2c17-11eb-3c0b-b79f49a927cb
@bind t_c Slider(1/alpha_t:0.1:4, default = 4, show_value=true)

# ╔═╡ 523e7660-2c16-11eb-10fd-0db54c6579f9
ϕ(x,t) = begin 
	r = sqrt((alpha_x*(x-x_c))^2 + (alpha_t*(t-t_c))^2) 
	return exp(-1/(1-r^2))*(r<1)
end

# ╔═╡ fad40bd0-2c18-11eb-38b5-81e9449bb96c
begin
	errors = Array{Float64,1}(undef, length(Ks))
	for (idx,K) in enumerate(Ks)
		I1 = integral(
			T, del_t(K), K, [t->emp_dens(K)(x,t) for x in 0:M-1], emp_dens(K), ϕ
		)
		I2 = integral(T, del_t(K), K, queue_fun, density_fun, ϕ)
		errors[idx] = abs(I1-I2)
	end
	plot(Ks, errors)
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
# ╠═c79f9e30-2855-11eb-2e82-2b9f4b94caf8
# ╟─9ccba9a0-285b-11eb-36da-b9b2eaad4a90
# ╠═a818baa0-285b-11eb-1766-99959ed836c0
# ╟─404d9450-2c1a-11eb-1837-4512f56f93ed
# ╠═e6413a10-2bed-11eb-00f6-5d88ebe4a8be
# ╠═26fc5480-2c0d-11eb-1aa8-f97eb83da977
# ╠═98d3d782-2c09-11eb-0875-79be1af36267
# ╠═2c206df0-2c0a-11eb-1a9f-a10b0b8ee4c7
# ╟─05f9e390-2c0a-11eb-1ce5-f1d4489d9273
# ╠═956e8b70-2c14-11eb-0268-3b5827f67b7f
# ╠═b7703960-2c16-11eb-3e0f-f3f194aa6f6a
# ╠═d4a9e9e0-2c16-11eb-1f56-8159f918e644
# ╠═8d024be0-2c17-11eb-1e1f-41b299015af9
# ╠═f73e3050-2c17-11eb-3c0b-b79f49a927cb
# ╠═523e7660-2c16-11eb-10fd-0db54c6579f9
# ╠═fad40bd0-2c18-11eb-38b5-81e9449bb96c
