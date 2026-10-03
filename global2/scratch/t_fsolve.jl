include("../geom.jl"); include("../rules.jl"); include("../fsolve.jl")
for n in (6, 8, 10)
    t = @elapsed s = fsolve([0.0,0,0,0], n)
    L = lg_bound(s.A0, s.A1, s.A2, 5π^2)
    println("square n=$n: th=", s.theta[1:3], " 2pi^2=", 2π^2, " 5pi^2=", 5π^2, " LG(beta=5pi^2)-2pi^2=", L-2π^2, " eps2=", s.eps2, " t=", t)
end
m = [0.25, 0.4330127]; nrm = [-0.8660254, 0.5]
for d in (0.03, 0.1, 0.2), (n,k) in ((10,10),(12,14),(10,14),(10,16))
    V = [0 0; 1 0; 0.5 0.8660254; (m+d*nrm)']
    s = fsolve_V(Float64.(V), n; kpsi=k)
    beta = 70.0
    KT = s.rho0 - s.eps2/(beta - s.rho0)
    L = lg_bound(s.A0, s.A1, s.A2, beta); L0 = lg_bound(s.A0, s.A1, s.A2_nopsi, beta)
    println("d=$d n=$n k=$k rho=", round(s.rho0,digits=4), " KT=", round(KT,digits=3), " LGnopsi=", round(L0,digits=3), " LG=", round(L,digits=4), " th2=", round(s.theta[2],digits=3))
end
