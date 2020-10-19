### A Pluto.jl notebook ###
# v0.12.4

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

# ╔═╡ f72e0a00-1208-11eb-3b00-337530f05d71
using DifferentialEquations

# ╔═╡ 552e37b0-1209-11eb-2fff-4b0cc5806afd
begin
	using Plots
	using PlutoUI
end

# ╔═╡ acc97ba0-1219-11eb-26cf-3999a572327e
using Zygote

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

# ╔═╡ 25efa4c0-1209-11eb-07b4-89eb9573557c
function follow_leader_solve(N, ρ_L, ρ_R)
    @assert iseven(N)

    velocity(ρ)=(1-ρ)
    length = 1/N

    x_0 = Array{Float64}(undef, N)
    for idx in 1:(Int32(N/2) -1)
        x_0[idx]=-1/(2*ρ_L) + (idx-1)/(N*ρ_L)
    end
    for idx in Int32(N/2):N
        x_0[idx] = idx/(N*ρ_R) - 1/(2*ρ_R)
    end

    return follow_leader_solve(velocity, length, N, ρ_R, x_0, (0.0,1.0))
end

# ╔═╡ 0cc137c0-1213-11eb-23ca-5b113e5aafb1
md"### N"

# ╔═╡ ebb1ee30-1212-11eb-0a1f-a11e6b2c4c38
@bind N Slider(10:2:200, default=100, show_value=true)

# ╔═╡ 5a2ed7b2-1213-11eb-00c8-d18c59de2236
md"### $\rho_L$"

# ╔═╡ da5a1f30-1213-11eb-15d0-c9d0aef75e19
@bind ρ_L Slider(0:0.05:1, default = 0.2, show_value=true)

# ╔═╡ 94b0b7a0-1213-11eb-0679-9786e2282683
md"### $\rho_R$"

# ╔═╡ a159fd3e-1213-11eb-1955-59d0a8180e96
@bind ρ_R Slider(0:0.05:1, default = 0.8, show_value=true)

# ╔═╡ 2ec1f1c0-1209-11eb-1cc2-adcc97ec8d9b
begin
    solution = follow_leader_solve(N, ρ_L, ρ_R)
    plot(solution, legend=false)
end

# ╔═╡ 92140930-1212-11eb-2d09-e5985aa73786
md"# Exercise 2"

# ╔═╡ 9ab249b0-1219-11eb-3a5b-917e0a0f82f5
import Pkg; Pkg.add("Zygote")

# ╔═╡ 66c49bb0-121b-11eb-36b4-a925d39c30f3
f(x) = x^2

# ╔═╡ 6dfa28a0-121b-11eb-18c0-9fa516f82579
f'(2)

# ╔═╡ b29dbd20-1214-11eb-2317-175b6d323e6b
""" assume flux concave """
function riemannProblem(density_left, density_right, flux, dflux, inv_dflux)
	if density_left < density_right
		s = (flux(density_right) - flux(density_left))/(density_right - density_left)
		density(x,t) = begin
			if x < s*t
				return density_left
			else
				return density_right
			end
		end
	else
		density(x,t) = begin
			if x <= dflux(density_left) *t
				return density_left
			elseif x >= dflux(density_right) * t
				return density_right
			else
				return inv_dflux(x/t)
			end
		end
	end
	return density
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
	dflux(ρ) = 1-2*ρ
	inv_dflux(x) = (1-x)/2
	exact_solution = riemannProblem(density_left, density_right, flux, dflux, inv_dflux)
end

# ╔═╡ 6b290da2-1223-11eb-1fbf-73eb07577b0d
@bind steppower Slider(1:5)

# ╔═╡ b0573c9e-1221-11eb-1f70-7dc5bdb9b80c
begin
	stepsize = 10^steppower
	times = 0:1/stepsize:1.0
	"times = $(times)"
end

# ╔═╡ e1a41620-1221-11eb-0d09-773c79dd9728
md"### Steps: $(stepsize)"

# ╔═╡ 47582120-1225-11eb-32f9-ab9eebb65767
plot(times, [exact_solution(0.5,t) for t in times])

# ╔═╡ a04023a0-1225-11eb-06ce-554f90454e1d
2.^1:10

# ╔═╡ 9cc1e970-1211-11eb-1833-81bd4c962dee
md"# Exercise 3"

# ╔═╡ 8d4929c0-1209-11eb-0751-7ba877da1f1f
HTML("<style> main { max-width:1000px; } </style> ")

# ╔═╡ Cell order:
# ╠═f72e0a00-1208-11eb-3b00-337530f05d71
# ╟─8d9d1520-1214-11eb-33f0-535efa40ff15
# ╠═229a5f40-1209-11eb-1f10-85c2da29fbd2
# ╠═25efa4c0-1209-11eb-07b4-89eb9573557c
# ╠═552e37b0-1209-11eb-2fff-4b0cc5806afd
# ╟─0cc137c0-1213-11eb-23ca-5b113e5aafb1
# ╟─ebb1ee30-1212-11eb-0a1f-a11e6b2c4c38
# ╟─5a2ed7b2-1213-11eb-00c8-d18c59de2236
# ╟─da5a1f30-1213-11eb-15d0-c9d0aef75e19
# ╟─94b0b7a0-1213-11eb-0679-9786e2282683
# ╟─a159fd3e-1213-11eb-1955-59d0a8180e96
# ╠═2ec1f1c0-1209-11eb-1cc2-adcc97ec8d9b
# ╟─92140930-1212-11eb-2d09-e5985aa73786
# ╠═9ab249b0-1219-11eb-3a5b-917e0a0f82f5
# ╠═acc97ba0-1219-11eb-26cf-3999a572327e
# ╠═66c49bb0-121b-11eb-36b4-a925d39c30f3
# ╠═6dfa28a0-121b-11eb-18c0-9fa516f82579
# ╠═b29dbd20-1214-11eb-2317-175b6d323e6b
# ╟─7559cea0-1222-11eb-26ea-a9621376ac85
# ╟─865ecd40-1222-11eb-0cdf-e109bea249c0
# ╟─e5587620-1222-11eb-03f3-f5aba45424c9
# ╟─ed920590-1222-11eb-3667-e1d9354db5b3
# ╠═9c554720-1220-11eb-39d1-6103b1cd7a3d
# ╟─e1a41620-1221-11eb-0d09-773c79dd9728
# ╟─6b290da2-1223-11eb-1fbf-73eb07577b0d
# ╟─b0573c9e-1221-11eb-1f70-7dc5bdb9b80c
# ╠═47582120-1225-11eb-32f9-ab9eebb65767
# ╠═a04023a0-1225-11eb-06ce-554f90454e1d
# ╟─9cc1e970-1211-11eb-1833-81bd4c962dee
# ╠═8d4929c0-1209-11eb-0751-7ba877da1f1f
