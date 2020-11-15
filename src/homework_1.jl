### A Pluto.jl notebook ###
# v0.12.7

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

# ╔═╡ 552e37b0-1209-11eb-2fff-4b0cc5806afd
begin
	using Plots
	using PlutoUI
end

# ╔═╡ f72e0a00-1208-11eb-3b00-337530f05d71
import DifferentialEquations: ODEProblem, solve

# ╔═╡ 8d9d1520-1214-11eb-33f0-535efa40ff15
md"# Exercise 1"

# ╔═╡ 229a5f40-1209-11eb-1f10-85c2da29fbd2
function follow_leader_solve(velocity, length , number, leader_density, x_0, tspan)
    function ode!(dx, x, p ,t)
        for idx in 1:(number-1)
            dx[idx] = velocity(length/(x[idx+1]-x[idx]))
        end
        dx[number] = velocity(leader_density)
    end
    problem = ODEProblem(ode!, x_0, tspan)
    return solve(problem, alg_hints=[:stiff])
    # return solve(problem, Rosenbrock23() #= ode23s equivalent =#)
end

# ╔═╡ 85444f00-1242-11eb-1f40-799386a407f1
function get_x_0(N, ρ_L, ρ_R)
	x_0 = Array{Float64}(undef, N)
    for idx in 1:(Int32(N/2) -1)
        x_0[idx]=-1/(2*ρ_L) + (idx-1)/(N*ρ_L)
    end
    for idx in Int32(N/2):N
        x_0[idx] = idx/(N*ρ_R) - 1/(2*ρ_R)
    end
	return x_0
end

# ╔═╡ 25efa4c0-1209-11eb-07b4-89eb9573557c
function follow_leader_solve(N, ρ_L, ρ_R, T)
    @assert iseven(N)

    velocity(ρ)=(1-ρ)
    length = 1/N

    x_0 = get_x_0(N, ρ_L, ρ_R)

    return follow_leader_solve(velocity, length, N, ρ_R, x_0, (0.0,T))
end

# ╔═╡ d9f8fc80-124c-11eb-37e1-df454b1614c9
md"### T"

# ╔═╡ e1ba0090-124c-11eb-1d0b-f3f6247c915f
@bind T_ex1 Slider(0.5:0.5:10, default=1, show_value=true)

# ╔═╡ 0cc137c0-1213-11eb-23ca-5b113e5aafb1
md"### N"

# ╔═╡ ebb1ee30-1212-11eb-0a1f-a11e6b2c4c38
@bind N Slider(10:2:200, default=100, show_value=true)

# ╔═╡ 5a2ed7b2-1213-11eb-00c8-d18c59de2236
md"### $\rho_L$"

# ╔═╡ da5a1f30-1213-11eb-15d0-c9d0aef75e19
@bind ρ_L Slider(0.05:0.05:1, default = 0.2, show_value=true)

# ╔═╡ 94b0b7a0-1213-11eb-0679-9786e2282683
md"### $\rho_R$"

# ╔═╡ a159fd3e-1213-11eb-1955-59d0a8180e96
@bind ρ_R Slider(0.05:0.05:1, default = 0.8, show_value=true)

# ╔═╡ d3c9a17e-1237-11eb-0629-898fb473e2d7
begin
	micro_sol = follow_leader_solve(N,ρ_L, ρ_R, T_ex1)
	p1 = plot()
	for idx in 1:length(micro_sol.u[1])
		plot!(p1, [x[idx] for x in micro_sol.u], micro_sol.t, legend=false)
	end
	p1
end

# ╔═╡ 92140930-1212-11eb-2d09-e5985aa73786
md"# Exercise 2"

# ╔═╡ 9ab249b0-1219-11eb-3a5b-917e0a0f82f5
import Pkg; Pkg.add("Zygote")

# ╔═╡ acc97ba0-1219-11eb-26cf-3999a572327e
import Zygote

# ╔═╡ b29dbd20-1214-11eb-2317-175b6d323e6b
""" assume flux concave """
function riemannProblem(density_left, density_right, flux, dflux, inv_dflux)
	if density_left < density_right
		s = (flux(density_right) - flux(density_left))/(density_right - density_left)
		return (x,t) -> begin
			if x < s*t
				return density_left
			else
				return density_right
			end
		end
	else
		return (x,t) -> begin
			if x <= dflux(density_left) *t
				return density_left
			elseif x >= dflux(density_right) * t
				return density_right
			else
				return inv_dflux(x/t)
			end
		end
	end
end

# ╔═╡ 7559cea0-1222-11eb-26ea-a9621376ac85
md"### $\rho_L$"

# ╔═╡ 865ecd40-1222-11eb-0cdf-e109bea249c0
@bind density_left Slider(0.05:0.05:1.0, default=0.2, show_value=true)

# ╔═╡ e5587620-1222-11eb-03f3-f5aba45424c9
md"### $\rho_R$"

# ╔═╡ ed920590-1222-11eb-3667-e1d9354db5b3
@bind density_right Slider(0.05:0.05:1.0, default = 0.8, show_value=true)

# ╔═╡ 9c554720-1220-11eb-39d1-6103b1cd7a3d
begin
	flux(ρ) = ρ*(1-ρ)
	inv_dflux(x) = (1-x)/2
	exact_solution = riemannProblem(density_left, density_right, flux, flux', inv_dflux)
end

# ╔═╡ 6b290da2-1223-11eb-1fbf-73eb07577b0d
@bind steppower Slider(1:5, default=2, show_value=true)

# ╔═╡ e1a41620-1221-11eb-0d09-773c79dd9728
begin
	stepsize = 10^steppower
	steps = -1.0:1/stepsize:1.0
	md"### Interpolation steps: $(2*stepsize)"
end

# ╔═╡ 49262c90-1234-11eb-1eab-d19a2c86da9e
@bind T Slider(1:0.5:5, default = 2, show_value=true)

# ╔═╡ 26d9a810-1234-11eb-09ee-8fbab86aba3d
md"until T: $(T)" 

# ╔═╡ 74208720-1229-11eb-3d0e-9f4373534945
@bind timesteps Slider(1:10, default=5, show_value=true)

# ╔═╡ 2790c410-1229-11eb-3ab9-37512e7dc4c3
begin
	times = 0:T/timesteps:T
	md"### Times: $(times)"
end

# ╔═╡ 9d0daa8e-1234-11eb-0739-2d340e62799f
md"Steps: $(timesteps)"

# ╔═╡ 47582120-1225-11eb-32f9-ab9eebb65767
begin
	p = plot()
	for t in times
		plot!(p, steps, [exact_solution(x,t) for x in steps], label = "t=$(t)", linetype=:steppost)
	end
	p
end

# ╔═╡ 15ec5bb0-122a-11eb-03ca-55cb5d64ae14
@bind single_time Slider(0:0.05:5, default = 0, show_value=true)

# ╔═╡ a04023a0-1225-11eb-06ce-554f90454e1d
plot(steps, [exact_solution(x,single_time) for x in steps], label="ρ(x,$(single_time))", linetype=:steppost, ylims=Tuple(sort([density_left, density_right])))

# ╔═╡ 9cc1e970-1211-11eb-1833-81bd4c962dee
md"# Exercise 3"

# ╔═╡ df746b70-124e-11eb-37d6-4148a9899621
Pkg.add("Polynomials")

# ╔═╡ e7081df0-124e-11eb-158f-fb7b47dfd45b
import Polynomials

# ╔═╡ 62594092-1242-11eb-2fc1-8b86a7c8beb9
function local_density(micro_solution, L)
	N = length(micro_solution[1])
	return [
		L/(x[idx+1]-x[idx]) 
		for idx in 1:N-1, x in micro_solution
	]
end

# ╔═╡ 979c5320-1245-11eb-13b0-7b09cfa70c0b
function l1_error(micro_solution, times, L, exact_solution)
	loc_density = local_density(micro_solution, L)
	
	timsteps = length(times)
	err = zeros(timesteps)
	for t in 1:timesteps
		x_t = micro_solution[t]
		time = times[t]
		for idx in 1:length(x_t)-1
			err[t] += abs(loc_density[idx,t] - exact_solution(x_t[idx],time))*(x_t[idx+1]-x_t[idx])
		end
	end
	return err
end

# ╔═╡ 70ced9d0-16b8-11eb-2cd3-295595c0fb33
md"### $\rho_L$"

# ╔═╡ 8f37bb80-16b8-11eb-1faf-4f8d7a5fc477
@bind rho_L Slider(0.05:0.05:1, default = 0.2, show_value=true)

# ╔═╡ 7dbb69fe-16b8-11eb-3b02-63b46e8bce9b
md"### $\rho_R$"

# ╔═╡ 96790840-16b8-11eb-3fbc-71c02308bf75
@bind rho_R Slider(0.05:0.05:1, default = 0.8, show_value=true)

# ╔═╡ 57441200-124e-11eb-2707-67b080d99934
md"### T"

# ╔═╡ 3f027700-124d-11eb-2fb6-5f4409ac3c82
@bind T_ex3 Slider(0.5:0.5:10, default=2, show_value=true)

# ╔═╡ c88f5000-1249-11eb-1428-bb7ca494fd97
begin
	n_vector = [50,100,200,300,400,800]
	micro_solutions = [follow_leader_solve(n, rho_L, rho_R, T_ex3) for n in n_vector]
	err = [l1_error(m_sol.u, m_sol.t, 1/n, exact_solution) for (n, m_sol) in zip(n_vector, micro_solutions)]
	plot(log.(n_vector), [log(x[T]) for x in err], label="log(N)-log(err[T])")
end

# ╔═╡ 1fabc800-124f-11eb-39b1-2708effb11c3
md"### degree"

# ╔═╡ 36cc3100-124f-11eb-071c-a3af8e47c1a2
@bind deg Slider(1:length(n_vector)-1, default=1, show_value=true)

# ╔═╡ f2171480-124e-11eb-1850-4f30b39bc919
Polynomials.fit(log.(n_vector), [log(x[T]) for x in err], deg)

# ╔═╡ 8d4929c0-1209-11eb-0751-7ba877da1f1f
HTML("<style> main { max-width:1000px; } </style> ")

# ╔═╡ Cell order:
# ╠═f72e0a00-1208-11eb-3b00-337530f05d71
# ╟─8d9d1520-1214-11eb-33f0-535efa40ff15
# ╠═229a5f40-1209-11eb-1f10-85c2da29fbd2
# ╠═85444f00-1242-11eb-1f40-799386a407f1
# ╠═25efa4c0-1209-11eb-07b4-89eb9573557c
# ╠═552e37b0-1209-11eb-2fff-4b0cc5806afd
# ╟─d9f8fc80-124c-11eb-37e1-df454b1614c9
# ╟─e1ba0090-124c-11eb-1d0b-f3f6247c915f
# ╟─0cc137c0-1213-11eb-23ca-5b113e5aafb1
# ╟─ebb1ee30-1212-11eb-0a1f-a11e6b2c4c38
# ╟─5a2ed7b2-1213-11eb-00c8-d18c59de2236
# ╠═da5a1f30-1213-11eb-15d0-c9d0aef75e19
# ╟─94b0b7a0-1213-11eb-0679-9786e2282683
# ╠═a159fd3e-1213-11eb-1955-59d0a8180e96
# ╠═d3c9a17e-1237-11eb-0629-898fb473e2d7
# ╟─92140930-1212-11eb-2d09-e5985aa73786
# ╠═9ab249b0-1219-11eb-3a5b-917e0a0f82f5
# ╠═acc97ba0-1219-11eb-26cf-3999a572327e
# ╠═b29dbd20-1214-11eb-2317-175b6d323e6b
# ╟─7559cea0-1222-11eb-26ea-a9621376ac85
# ╠═865ecd40-1222-11eb-0cdf-e109bea249c0
# ╟─e5587620-1222-11eb-03f3-f5aba45424c9
# ╠═ed920590-1222-11eb-3667-e1d9354db5b3
# ╠═9c554720-1220-11eb-39d1-6103b1cd7a3d
# ╟─e1a41620-1221-11eb-0d09-773c79dd9728
# ╟─6b290da2-1223-11eb-1fbf-73eb07577b0d
# ╠═2790c410-1229-11eb-3ab9-37512e7dc4c3
# ╟─26d9a810-1234-11eb-09ee-8fbab86aba3d
# ╠═49262c90-1234-11eb-1eab-d19a2c86da9e
# ╟─9d0daa8e-1234-11eb-0739-2d340e62799f
# ╠═74208720-1229-11eb-3d0e-9f4373534945
# ╠═47582120-1225-11eb-32f9-ab9eebb65767
# ╠═15ec5bb0-122a-11eb-03ca-55cb5d64ae14
# ╠═a04023a0-1225-11eb-06ce-554f90454e1d
# ╟─9cc1e970-1211-11eb-1833-81bd4c962dee
# ╠═df746b70-124e-11eb-37d6-4148a9899621
# ╠═e7081df0-124e-11eb-158f-fb7b47dfd45b
# ╠═62594092-1242-11eb-2fc1-8b86a7c8beb9
# ╠═979c5320-1245-11eb-13b0-7b09cfa70c0b
# ╟─70ced9d0-16b8-11eb-2cd3-295595c0fb33
# ╟─8f37bb80-16b8-11eb-1faf-4f8d7a5fc477
# ╟─7dbb69fe-16b8-11eb-3b02-63b46e8bce9b
# ╟─96790840-16b8-11eb-3fbc-71c02308bf75
# ╟─57441200-124e-11eb-2707-67b080d99934
# ╟─3f027700-124d-11eb-2fb6-5f4409ac3c82
# ╠═c88f5000-1249-11eb-1428-bb7ca494fd97
# ╟─1fabc800-124f-11eb-39b1-2708effb11c3
# ╟─36cc3100-124f-11eb-071c-a3af8e47c1a2
# ╠═f2171480-124e-11eb-1850-4f30b39bc919
# ╟─8d4929c0-1209-11eb-0751-7ba877da1f1f
