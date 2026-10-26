using Test, Random, LinearAlgebra
include(joinpath(@__DIR__, "run_revision_sensitivity.jl"))
const R = RevisionSensitivity
const B = R.B

function brute_repair(flow, qubo)
    terms = [(i,j,qubo.quadratic[i,j]) for i in 1:36 for j in 1:36
             if qubo.quadratic[i,j] != 0]
    x = vcat(flow, zeros(Int,17))
    best = Inf
    for mask in 0:(1 << 17)-1
        for i in 1:17
            x[19+i] = (mask >> (i-1)) & 1
        end
        value = qubo.offset + dot(qubo.linear,x)
        for (i,j,c) in terms
            value += c*x[i]*x[j]
        end
        best = min(best,qubo.scale*value)
    end
    return best
end

@testset "Exact repair and penalty scaling" begin
    q,c,t = R.load_problem()
    rng = MersenneTwister(97191)
    for y in (B.flow_from_index(parse(Int,R.REFERENCE;base=2),19), B.random_flow(rng,19))
        nominal = B.repair_flow(y,q,c).exact_repaired_qubo_energy
        @test nominal ≈ brute_repair(y,q) atol=1e-8
        for alpha in (.01,.25,4.)
            scaled = B.QuboData(q.scale, alpha*q.offset,
                alpha*q.linear + (1-alpha)/q.scale*vcat(R.COST,zeros(17)),alpha*q.quadratic)
            predicted = dot(R.COST,y)+alpha*(nominal-dot(R.COST,y))
            @test predicted ≈ brute_repair(y,scaled) atol=1e-8
        end
    end
end

@testset "Walsh fit equals ordinary least squares" begin
    rng = MersenneTwister(97192)
    values = randn(rng,256)
    coefficients,masks = R.fit_full(values,8)
    X = R.design_matrix(1:256,masks)
    @test coefficients ≈ X\values atol=1e-12
    @test R.predict(coefficients,masks,8) ≈ X*coefficients atol=1e-12
    q,c,t = R.load_problem()
    energy,cost = R.enumerate_objective(q,c)
    a,m = R.fit_full(energy,19)
    prediction = R.predict(a,m,19)
    reduced = joinpath(R.ROOT,"ds_mfg_reduced_flow_objective")
    scalars = split(readlines(joinpath(reduced,"reduced_scalars.csv"))[2],',')
    scale,offset = parse.(Float64,scalars[2:3])
    L = [parse(Float64,split(s,',')[2]) for s in readlines(joinpath(reduced,"reduced_L_vector.csv"))[2:end]]
    Q = reduce(vcat,[permutedims(parse.(Float64,split(s,',')[2:end])) for s in
        readlines(joinpath(reduced,"reduced_Q_matrix.csv"))[2:end]])
    for bits in keys(t)
        mask = parse(Int,bits;base=2)
        y = B.flow_from_index(mask,19)
        @test prediction[mask+1] ≈ scale*(offset+dot(L,y)+dot(y,Q*y)) atol=1e-8
    end
    rows = R.penalty_audit(energy,cost,t)
    @test !first(rows).optimum_feasible
    @test only(filter(r->r.penalty_multiplier==1.,rows)).min_infeasible_minus_max_feasible > 0
end

@testset "Search budget and seed reproducibility" begin
    q,c,t = R.load_problem()
    for method in ("uniform","hill_climb","simulated_annealing")
        a = R.baseline_run(method,97193,1003,q,c)
        b = R.baseline_run(method,97193,1003,q,c)
        @test a.evaluations == 1003
        @test a.unique_evaluations <= a.evaluations
        @test a.best_repaired_energy == b.best_repaired_energy
        @test a.optimum_evaluations == b.optimum_evaluations
        @test a.restarts == b.restarts
    end
end
