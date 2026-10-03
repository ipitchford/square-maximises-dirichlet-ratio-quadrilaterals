include("../cover_t.jl")
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
function kpt(K, M, A2, beta, L0)
    th = eigvals(Symmetric(K), Symmetric(M)); L = beta*M[1,1] > K[1,1] ? (beta*M[1,1]-K[1,1])/(beta*A2-M[1,1]) : 0.0
    th[2]/max(L, L0), th
end
lines = filter(l -> occursin("pilot_pass", l), readlines(joinpath(@__DIR__, "..", "runs", "pilot_t1", "cell_1.jsonl")))
rng = MersenneTwister(2); W = TWorker()
for ln in lines[1:6]
    lo = parse.(Float64, split(match(r"\"lo\":\[([^\]]+)\]", ln).captures[1], ",")); hi = parse.(Float64, split(match(r"\"hi\":\[([^\]]+)\]", ln).captures[1], ","))
    n = parse(Int, match(r"\"n\":(\d+)", ln).captures[1])
    # parent box: double the box in its last split dimension (approximate: double widest relative)
    p0, rad = box_center_radius(lo, hi)
    for scale in (1.0, 2.0, 4.0)
        r = rad .* scale
        sb = scalar_bounds(p0, r); L0 = max(sb.L_rect, sb.L_hersch); beta = max(sb.beta_rect_ks, addlo(L0, sb.gapAC))
        line = @sprintf("p0=%s rad=%s x%.0f beta=%.2f:", string(round.(p0,digits=3)), string(round.(rad,digits=4)), scale, beta)
        for sh in (0, 1)
            s = tsolve(p0, n; shift = sh); R = transport_refs(Float64, s)
            b = transport_base!(W.tpF, R, p0, r, W.wsF, sh)
            res = crit_eval!(W.cpF, W.tpF.M, W.tpF.K, W.tpF.A2, W.wsF; mode = :LG, beta, L0)
            kc = res.ok ? kappa_star_t!(W.cpF, (; M = W.tpF.M, K = W.tpF.K, A2 = W.tpF.A2), W.wsF, :LG, beta, L0; steps = 10) : NaN
            kt = 0.0; ktb = 0.0; th2min = Inf
            for smp in 1:40
                t = smp <= 16 ? [((smp-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
                Kp, Mp, A2p = tpencil(R, p0 .+ r .* t, sh)
                k, th = kpt(Kp, Mp, A2p, beta, L0); kt = max(kt, k); th2min = min(th2min, th[2])
            end
            for smp in 1:40   # same with an (optimistic) beta = 0.9 * min theta2 over the box
                t = smp <= 16 ? [((smp-1)>>k)&1 == 1 ? 1.0 : -1.0 for k in 0:3] : 2 .* rand(rng, 4) .- 1
                Kp, Mp, A2p = tpencil(R, p0 .+ r .* t, sh)
                k, th = kpt(Kp, Mp, A2p, 0.9*th2min, L0); ktb = max(ktb, k)
            end
            line *= @sprintf("  sh%d: cert %s true %.3f true(beta=0.9th2) %.3f", sh, isnan(kc) ? "FAIL " : @sprintf("%.3f", kc), kt, ktb)
        end
        println(line)
    end
end
