# tests.jl — consolidated self-tests of the global2 kernel (run: julia tests.jl; ≈ 1–2 min)
include("cover.jl")
using Random, Test

@testset "global2 kernel" begin
    # 1. square: RR eigenvalues and LG bound
    s = fsolve([0.0, 0, 0, 0], 10)
    @test abs(s.theta[1] - 2π^2) < 1e-11 && abs(s.theta[2] - 5π^2) < 1e-9
    @test abs(lg_bound(s.A0, s.A1, s.A2, 5π^2) - 2π^2) < 1e-10
    # 2. point box reproduces the solver (M = I, K = diag θ, A2 via the identity)
    p0 = [0.05, 0.2, 0.1, -0.05]
    sd = fsolve_d(p0, 8)
    aw = newasm(Float64)
    b0 = assemble_box!(aw, p0, zeros(4), sd)
    @test maximum(abs, [x.c[1] for x in b0.M] .- [1, 0, 0, 1, 0, 1]) < 1e-12
    @test maximum(abs, [x.c[1] for x in b0.K] .- [sd.theta[1], 0, 0, sd.theta[2], 0, sd.theta[3]]) < 1e-9
    @test abs(b0.A2.c[1] - sd.A2) / sd.A2 < 1e-10
    # 3. Taylor model encloses exact values at random points of a box (Float64 check of the formulas)
    rad = [0.01, 0.02, 0.015, 0.01]
    bB = assemble_box!(aw, p0, rad, sd)
    TMs = vcat(deepcopy(bB.M), deepcopy(bB.K), [deepcopy(bB.A2)])
    aw2 = newasm(Float64); rng = MersenneTwister(1); worst = 0.0
    for _ = 1:15
        t = 2 .* rand(rng, 4) .- 1
        bp = assemble_box!(aw2, p0 .+ rad .* t, zeros(4), sd)
        vals = vcat([x.c[1] for x in bp.M], [x.c[1] for x in bp.K], [bp.A2.c[1]])
        for q = 1:13
            worst = max(worst, abs(tm_eval(TMs[q], t) - vals[q]) / remval(TMs[q]))
        end
    end
    @test worst < 1.0
    # 4. Arb Taylor models contain the Float64 ones (coefficients within radius + rounding slack, equal remainders)
    awA = newasm(Arb)
    bA = assemble_box!(awA, p0, rad, sd)
    for (fA, fF) in zip(vcat(bA.M, bA.K, [bA.A2]), TMs)
        scale = sum(abs, fF.c)        # Float64 rounding noise is relative to the TM's magnitude, not to each coefficient
        for k = 1:TAB[].nc
            @test abs(Float64(Arblib.midref(fA.c[k])) - fF.c[k]) <= Float64(Arblib.radref(fA.c[k])) + 1e-12 * scale
        end
        @test abs(remval(fA) - remval(fF)) <= 1e-6 * remval(fF) + 1e-300
    end
    # 5. criterion: certifies a moderate box rigorously, not an inflated one
    cpA = critpool(Arb); wsA = newwork(Arb)
    sb = scalar_bounds(p0, rad); L0 = max(sb.L_rect, sb.L_hersch); beta = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
    r = crit_eval!(cpA, bA.M, bA.K, bA.A2, wsA; mode = :LG, beta, L0)
    @test r.ok
    rbad = crit_eval!(cpA, bA.M, bA.K, bA.A2, wsA; mode = :LG, beta, L0, kappa = 3 // 2)
    @test !rbad.ok          # ξ ≈ 2.3 here, so ξ < 3/2 must not be certified
    # 6. exclusions: the local cube, a non-convex box, a near-triangle box
    @test exclusion([0.0, 0.0, 0.0, -0.01], [0.01, 0.01, 0.01, 0.01]; cube = CUBE) === :E1
    @test exclusion([0.0, 2.5, 2.5, 0.0], [0.01, 2.51, 2.51, 0.01]; cube = CUBE) in (:E4, :E2, :E3)
    # 7. D4 invariance of the Float64 pipeline
    p = [0.16, 0.483, 0.432, 0.422]
    th = [fsolve(q, 10).theta[2] for q in (p, rho_rot(p), sigma_h(p), sigma_v(p))]
    @test maximum(th) - minimum(th) < 1e-8 * th[1]
    # 8. box containment is exact
    p0c, rc = box_center_radius([0.1, 0.2, 0.3, -0.4], [0.3, 0.7, 0.30000000000000004, 0.1])
    @test all(Rational{BigInt}(p0c[i]) - Rational{BigInt}(rc[i]) <= Rational{BigInt}([0.1, 0.2, 0.3, -0.4][i]) for i = 1:4)
end
