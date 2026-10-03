# Prototype: transported 1-D lower bound λ1(Q_p) ≥ c(box)·μ1D(Q_p0) on E3T leaf boxes of the transport cover.
include("../replay_t.jl"); include("../bound1d.jl")
using Printf, Random
function transport_factor(p0, rad, sh)
    perm = vperm(sh)
    VA = [let c = Arb(V0[perm[j], k]; prec = PREC)
              for i = 1:4; c += Arb(DV[perm[j], k, i]; prec = PREC) * Arb(p0[i]; prec = PREC); end
              c
          end for j = 1:4, k = 1:2]
    amin = Inf; bmax = 0.0
    for (i0, i1, i2, _) in TRIS
        M0 = [VA[i1, 1]-VA[i0, 1] VA[i2, 1]-VA[i0, 1]; VA[i1, 2]-VA[i0, 2] VA[i2, 2]-VA[i0, 2]]
        d0 = M0[1, 1] * M0[2, 2] - M0[1, 2] * M0[2, 1]
        Arblib.contains_zero(d0) && return 0.0           # orientation may be negative (TRIS sign)
        Mi = [M0[2, 2]/d0 -M0[1, 2]/d0; -M0[2, 1]/d0 M0[1, 1]/d0]
        eab = [Arb(0; prec = PREC) for _ = 1:2, _ = 1:2]; trc = Arb(0; prec = PREC)
        for i = 1:4
            dM = [Arb(DV[perm[i1], 1, i] - DV[perm[i0], 1, i]; prec = PREC) Arb(DV[perm[i2], 1, i] - DV[perm[i0], 1, i]; prec = PREC);
                  Arb(DV[perm[i1], 2, i] - DV[perm[i0], 2, i]; prec = PREC) Arb(DV[perm[i2], 2, i] - DV[perm[i0], 2, i]; prec = PREC)]
            C = dM * Mi
            r = Arb(rad[i]; prec = PREC)
            for a = 1:2, b = 1:2; eab[a, b] += abs(C[a, b]) * r; end
            trc += abs(C[1, 1] + C[2, 2]) * r
        end
        eF = sqrt(eab[1, 1]^2 + eab[1, 2]^2 + eab[2, 1]^2 + eab[2, 2]^2)
        dE = eab[1, 1] * eab[2, 2] + eab[1, 2] * eab[2, 1]
        dlo = 1 - trc - dE; dhi = 1 + trc + dE
        (Arblib.ispositive(1 - eF) && Arblib.ispositive(dlo)) || return 0.0
        a1 = (1 - eF) / (1 + eF); a2 = dlo / (1 + eF)^2
        al = max(Float64(Arblib.lbound(a1), RoundDown), Float64(Arblib.lbound(a2), RoundDown))
        be = min(Float64(Arblib.ubound(dhi), RoundUp), Float64(Arblib.ubound((1 + eF)^2), RoundUp))
        amin = min(amin, al); bmax = max(bmax, be)
    end
    return Float64(Arblib.lbound(Arb(amin; prec = PREC) / Arb(bmax; prec = PREC)), RoundDown)
end
logf = ARGS[1]; K = parse(Int, ARGS[2])
recs = [r for r in (parse_rec(ln) for ln in eachline(logf)) if get(r, "status", "") == "E3T"]
rng = MersenneTwister(7)
sel = recs[randperm(rng, length(recs))[1:min(K, end)]]
@printf("%-5s %-8s %-6s %-6s | per r: c, T1D/θ1, monoκ, β0/θ1\n", "w/D", "θ1", "L1D/θ", "ξRR")
for r in sel
    p0, _ = box_center_radius(r["lo"], r["hi"])
    minimum(cvals(p0)) > 0 || continue
    wd = width_diam(verts(p0))
    th = tsolve1(p0, 8).theta
    L = lambda1_1d(p0, zeros(4); K = 48)
    line = @sprintf("%.3f %8.2f %.3f %.3f |", wd[1] / wd[2], th[1], L / th[1], th[2] / th[1])
    for rr in (0.005, 0.01, 0.02, 0.04)
        rad = fill(rr, 4)
        c = max(transport_factor(p0, rad, 0), transport_factor(p0, rad, 1))
        sb = scalar_bounds(p0, rad); L0 = max(L * c, sb.L_rect, sb.L_hersch)
        line *= @sprintf(" r=%.3f %.3f %.3f %.2f %.3f |", rr, c, L * c / th[1], th[2] / L0, (L0 + sb.gapAC) / th[1])
    end
    println(line)
end
