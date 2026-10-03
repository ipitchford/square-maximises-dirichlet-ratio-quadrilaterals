using Arblib
setprecision(Arb, 128)
x = Arb(0.3); m = Mag(0)
for mm in methods(Arblib.get!); println(mm); end
println(methods(Arblib.get_mag_lower!))
# time the right one
Arblib.get!(m, x)
t = @elapsed for i in 1:100000; Arblib.get!(m, x); end
println("get! ns: ", t/1e5*1e9)
f(m,x,N) = (for i in 1:N; Arblib.get!(m, x); end)
f(m,x,10); t=@elapsed f(m,x,100000); println("in func ns: ", t/1e5*1e9)
println(m)
