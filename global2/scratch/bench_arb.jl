using Arblib, LinearAlgebra
setprecision(Arb, 128)
println(pkgversion(Arblib))
x = [Arb(rand(); prec=128) for _ in 1:16]; y = [Arb(rand(); prec=128) for _ in 1:16]; z = [Arb(0; prec=128) for _ in 1:16]
function f!(z, x, y, N)
    for it in 1:N
        for k in 1:15
            Arblib.addmul!(z[k], x[k], y[16-k])
        end
    end
end
f!(z,x,y,10)
t = @elapsed f!(z, x, y, 100000)
println("Arb addmul! ns/op: ", t/(1.5e6)*1e9)
m = Mag(0); m2 = Mag(1.5); m3 = Mag(2.5)
function g!(m, m2, m3, N)
    for it in 1:N
        Arblib.addmul!(m, m2, m3)
    end
end
g!(m,m2,m3,10)
t = @elapsed g!(m, m2, m3, 1_000_000)
println("Mag addmul ns/op: ", t/1e6*1e9)
xf = rand(16); yf = rand(16); zf = zeros(16)
function h!(z, x, y, N)
    for it in 1:N
        @inbounds for k in 1:15
            z[k] += x[k]*y[16-k]
        end
    end
end
h!(zf,xf,yf,10)
t = @elapsed h!(zf, xf, yf, 1_000_000)
println("Float64 ns/op: ", t/1.5e7*1e9)
a = Arb(0.3); 
t = @elapsed for i in 1:100000; Arblib.get!(m, x[1]); end
println("get_mag ns: ", t/1e5*1e9)
println(Threads.nthreads())
