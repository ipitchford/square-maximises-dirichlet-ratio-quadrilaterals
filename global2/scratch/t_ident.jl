include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl"); include("../tm.jl"); include("../assemble.jl"); include("../criterion.jl")
set_order!(3)
for (p0, n) in (([0.05, 0.2, 0.1, -0.05], 8), ([0.2, 0.5, 0.1, -0.3], 12), ([0.5,0.3,0.2,0.3], 10))
    sol = fsolve(p0, n); aw = newasm(Float64)
    b = assemble_box!(aw, p0, zeros(4), sol)
    println("p0=$p0: A2 fsolve=", sol.A2, "  A2 identity=", b.A2.c[1], " rel diff=", (b.A2.c[1]-sol.A2)/sol.A2)
end
