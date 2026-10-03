include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl")
mm = [0.25, 0.4330127]; nrm = [-0.8660254, 0.5]
fem = Dict(0.03=>49.9826, 0.1=>44.8164, 0.2=>39.3150)
for dd in (0.03, 0.1, 0.2), n in (8, 10, 12, 14)
    V = Float64.([0 0; 1 0; 0.5 0.8660254; (mm+dd*nrm)'])
    t = @elapsed s = fsolve_V(V, n)
    L = lg_bound(s.A0, s.A1, s.A2, 70.0)
    println("d=$dd n=$n kpsi=$(s.kpsi): rho=", round(s.rho0,digits=4), " LG=", round(L,digits=4), " FEM=", fem[dd], " rel.gap=", round((fem[dd]-L)/fem[dd],sigdigits=3), "  th2=", round(s.theta[2],digits=4), " t=", round(t,digits=3))
end
for n in (6,8,10)
    s = fsolve([0.0,0,0,0], n); println("square n=$n LG-2pi^2 = ", lg_bound(s.A0,s.A1,s.A2,5π^2)-2π^2)
end
