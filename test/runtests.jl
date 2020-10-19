using scalarConservationLaws
using Test
using Plots


@testset "Exercise 1" begin
    N = 20
    ρ_L = 0.2
    ρ_R = 0.8

    solution = scalarConservationLaws.follow_leader_solve(N, ρ_L, ρ_R)
    plot(solution, legend=:outertopright)
end