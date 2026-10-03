include("../cover.jl"); include("../transport.jl")
using Printf, Random, LinearAlgebra
# direct Float64 transported pencil at p (for checking the Taylor models)
function tpencil(R::TRefs{Float64}, p0, p, shift)
    V = verts(p)[vperm(shift), :]; K = zeros(3,3); M = zeros(3,3); A2 = 0.0
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        Mp = [V[i1,1]-V[i0,1] V[i2,1]-V[i0,1]; V[i1,2]-V[i0,2] V[i2,2]-V[i0,2]]
        A = Mp * R.Minv0[t]; d = det(A); adj = [A[2,2] -A[1,2]; -A[2,1] A[1,1]]; Ga = adj*adj'; GA = A'*A
        for i in 1:3, j in 1:3; K[i,j] += sum(Ga[a,b]*R.I[t][a,b,i,j] for a in 1:2, b in 1:2)/d; end
        M .+= d .* R.M[t]
        A2 += d*sum(GA .* R.TT[t]) + 2*sum(GA .* R.TC[t]) + sum(GA .* R.CC[t])/d
    end
    K, M, A2
end
function kappa_pt(K, M, A2, beta, L0)
    th = eigvals(Symmetric(K), Symmetric(M)); L = beta*M[1,1] > K[1,1] ? (beta*M[1,1]-K[1,1])/(beta*A2-M[1,1]) : 0.0
    th[2]/max(L, L0)
end
function kappa_box(tp, cp, R, p0, rad, ws, beta, L0, shift; mode=:LG)
    b = transport_base!(tp, R, p0, rad, ws, shift)
    b.ok || return NaN
    r = crit_eval!(cp, tp.M, tp.K, tp.A2, ws; mode, beta, L0, kappa = 5//2)
    r.ok || return NaN
    lo, hi = 1.0, 2.5
    for _ in 1:12
        mid = (lo+hi)/2
        rr = crit_eval!(cp, tp.M, tp.K, tp.A2, ws; mode, beta, L0, kappa = rationalize(mid; tol=1e-9))
        rr.ok ? (hi = mid) : (lo = mid)
    end
    hi
end
set_order!(3)
tp = tpool(Float64); cp = critpool(Float64); ws = newwork(Float64)
rng = MersenneTwister(4)
pts = Dict("generic" => [0.3,0.2,0.4,0.1], "nearflat_a" => [0.171, 0.374, 0.092, -0.438], "nearflat_b" => [0.262, 0.31, 0.297, -0.291],
           "obtuse148" => [0.2398, 0.1997, 0.1621, -0.2382], "nearflat_c" => [0.211, 0.446, 0.262, -0.162])
for (name, p0) in sort(collect(pts))
    for n in (8, 10)
        t0 = time(); s = tsolve(p0, n); R = transport_refs(Float64, s); tsl = time() - t0
        K, M, A2 = tpoint(R)
        sb = scalar_bounds(p0, fill(1e-12, 4)); L0 = max(sb.L_rect, sb.L_hersch); beta = max(sb.beta_rect_ks, L0 + sb.gapAC)
        kt = kappa_pt(K, M, A2, beta, L0)
        f = fsolve(p0, n); kf = f.theta[2]/max(lg_bound(f.A0, f.A1, f.A2, beta), L0)
        line = @sprintf("%-10s n=%2d: point kappa transport %.6f vs fixed %.6f (A2 %.8f vs %.8f) [%.3fs] |", name, n, kt, kf, A2, f.A2, tsl)
        for r in (0.01, 0.02, 0.05)
            rad = fill(r, 4)
            sbb = scalar_bounds(p0, rad); L0b = max(sbb.L_rect, sbb.L_hersch); betab = max(sbb.beta_rect_ks, L0b + sbb.gapAC)
            kb = kappa_box(tp, cp, R, p0, rad, ws, betab, L0b, s.shift)
            # true sampled (transported pencil, same beta)
            kt2 = 0.0; maxerr = 0.0
            b = transport_base!(tp, R, p0, rad, ws, s.shift)
            for smp in 1:30
                t = smp <= 16 ? [((smp-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
                Kp, Mp, A2p = tpencil(R, p0, p0 .+ rad .* t, s.shift)
                kt2 = max(kt2, kappa_pt(Kp, Mp, A2p, betab, L0b))
                # TM check: D*K11 at t vs direct
                d1 = tm_eval(tp.d[1], t); d2 = tm_eval(tp.d[2], t)
                maxerr = max(maxerr, abs(tm_eval(tp.K[1], t) - d1*d2*Kp[1,1]) / max(remval(tp.K[1]), 1e-300))
            end
            line *= @sprintf(" r=%.2f: cert %s true %.4f (TMerr/rem %.2f)", r, isnan(kb) ? "FAIL  " : @sprintf("%.4f", kb), kt2, maxerr)
        end
        println(line, "  shift=", s.shift, " angles=", round.(angles_deg(verts(p0)), digits=1))
    end
end
