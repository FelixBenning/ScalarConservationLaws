module scalarConservationLaws

using DifferentialEquations

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

end # module
