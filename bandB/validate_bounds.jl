# validate_bounds.jl — compare the rigorous 1-D bounds (point boxes) with the 2-D FEM eigenvalues of validate_fem.py.
# Author: Anonymous.  usage: julia validate_bounds.jl [validation_fem.csv]  -> validation_bounds.log (stdout)
# Checks per shape: E ≤ λ1_FEM (the lower bound is below the FEM value) and R2 ≥ λ2_FEM (FEM λ2 is itself an upper
# approximation, so R2 ≥ λ2_FEM implies R2 ≥ λ2), with the relative margins printed.
include(joinpath(@__DIR__, "cert.jl"))
using Printf
f = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "validation_fem.csv")
lines = readlines(f)
hdr = split(lines[1], ",")
@printf("%-5s %-14s %-9s %-3s %-12s %-12s %-9s %-12s %-12s %-9s %-8s %-8s %s\n", "chart", "tag", "s", "LA", "E", "lam1_FEM",
        "l1/E-1", "R2", "lam2_FEM", "R2/l2-1", "xi_FEM", "R2/E", "cert")
nbad = 0; m1 = Inf; m2 = Inf; n = 0; worst = 0.0; ncert = 0; byc = Dict(1 => 0, 2 => 0)
for ln in lines[2:end]
    c = split(ln, ",")
    chart = parse(Int, c[1]); tag = c[2]
    p = ntuple(k -> (q = split(c[2+k], "/"); QQ(parse(BigInt, q[1]), parse(BigInt, q[2]))), 4)
    s = parse(Float64, c[7]); la = c[8]
    l1 = parse(Float64, c[9]); l2 = parse(Float64, c[10])
    B = Box(chart, p, p)
    μ1f, R2f, bd, F = float_eval(B)
    rb = rigorous_bounds(B, DEFAULT_PARAMS, bd, F)
    E = rb.E; R2 = up(rb.R2A)
    cert = certify(B).ok
    ok1 = E <= l1; ok2 = R2 >= l2
    global nbad += !(ok1 && ok2); global n += 1; byc[chart] += 1
    global m1 = min(m1, l1 / E - 1); global m2 = min(m2, R2 / l2 - 1); global worst = max(worst, R2 / E); global ncert += cert
    @printf("%-5d %-14s %-9.5f %-3s %-12.5f %-12.5f %-9.5f %-12.5f %-12.5f %-9.5f %-8.5f %-8.5f %s%s\n", chart, tag, s, la, E, l1,
            l1 / E - 1, R2, l2, R2 / l2 - 1, l2 / l1, R2 / E, cert, (ok1 && ok2) ? "" : "   <-- BRACKET VIOLATED")
end
@printf("\nSUMMARY: %d shapes (chart I: %d, chart II: %d); bracket violations: %d; min (λ1_FEM/E − 1) = %.4f; min (R2/λ2_FEM − 1) = %.4f; max R2/E = %.4f; point certificates %d/%d\n",
        n, byc[1], byc[2], nbad, m1, m2, worst, ncert, n)
