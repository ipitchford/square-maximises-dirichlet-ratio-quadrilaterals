include("../cover_t.jl"); include("../bounds_tm.jl")
using Printf, Random, LinearAlgebra
function tpencil(R::TRefs{Float64}, p, shift)
    V = verts(p)[vperm(shift), :]; K = zeros(3,3); M = zeros(3,3); A2 = 0.0
    for (t, (i0, i1, i2, sg)) in enumerate(TRIS)
        Mp = [V[i1,1]-V[i0,1] V[i2,1]-V[i0,1]; V[i1,2]-V[i0,2] V[i2,2]-V[i0,2]]
        A = Mp * R.Minv0[t]; d = det(A); adj = [A[2,2] -A[1,2]; -A[2,1] A[1,1]]; Ga = adj*adj'; GA = A'*A
        for i in 1:3, j in 1:3; K[i,j] += sum(Ga[a,b]*R.I[t][a,b,i,j] for a in 1:2, b in 1:2)/d; end
        M .+= d .* R.M[t]; A2 += d*sum(GA .* R.TT[t]) + 2*sum(GA .* R.TC[t]) + sum(GA .* R.CC[t])/d
    end
    K, M, A2
end
function kfrom(K, M, A2, beta, L0)
    th = eigvals(Symmetric(K), Symmetric(M)); L = beta*M[1,1] > K[1,1] ? (beta*M[1,1]-K[1,1])/(beta*A2-M[1,1]) : 0.0
    th[2]/max(L, L0)
end
W = TWorker(); rng = MersenneTwister(9)
for p0 in ([0.321, 0.606, 1.013, 0.668], [0.3,0.2,0.4,0.1], [0.171, 0.374, 0.092, -0.438], [0.5, 0.3, 0.2, 0.3], [0.343, 0.671, 1.061, 0.712])
    sh = choose_shift(p0); s = tsolve(p0, 8; shift = sh); R = transport_refs(Float64, s)
    line = @sprintf("p0=%-26s sh=%d", string(p0), sh)
    for r in (0.01, 0.02, 0.04)
        rad = fill(r, 4); bt = bound_tms(Float64, p0, rad, W.wsF)
        kfix = 0.0; kopt = 0.0
        for smp in 1:24
            t = smp <= 16 ? [((smp-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
            p = p0 .+ rad .* t
            minimum(cvals(p)) > 0 || continue
            beta = maximum(tm_eval(B, t) - remval(B) for (nm, B) in bt.betas; init = 0.0)
            L0 = maximum(tm_eval(L, t) - remval(L) for (nm, L) in bt.L0s; init = 0.0)
            Kp, Mp, A2p = tpencil(R, p, sh); kfix = max(kfix, kfrom(Kp, Mp, A2p, beta, L0))
            so = fsolve(p, 8)   # optimal test functions at p (full RR + LS flux at p)
            LG = beta*so.A1 > so.A0 ? lg_bound(so.A0, so.A1, so.A2, beta) : 0.0
            kopt = max(kopt, so.theta[2]/max(LG, L0))
        end
        line *= @sprintf("  r=%.2f: fixed-transport %.3f | optimal-at-p %.3f", r, kfix, kopt)
    end
    println(line)
end
