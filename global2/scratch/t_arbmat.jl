using Arblib
setprecision(Arb, 128)
A = [Arb(rand()) for i in 1:15, j in 1:400]; B = [Arb(rand()) for i in 1:400, j in 1:15]
AM = ArbMatrix(A; prec = 128); BM = ArbMatrix(B; prec = 128)
C = AM * BM
println(typeof(C), " ", size(C), " ", C[1,1])
t = @elapsed (for _ in 1:10; C = AM * BM; end)
println("15x400x15 arb_mat_mul: ", t/10*1e3, " ms")
Cm = Matrix{Arb}(undef, 15, 15); for i in 1:15, j in 1:15; Cm[i,j] = Arb(C[i,j]); end
println(Cm[2,3])
# elementwise scaling cost
t = @elapsed (for _ in 1:10; for j in 1:400, i in 1:15; Arblib.mul!(A[i,j], A[i,j], Arb(1.0)); end; end)
println("elementwise 6000 mul (alloc const): ", t/10*1e3, " ms")
