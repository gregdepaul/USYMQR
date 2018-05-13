include("tridiagonalization.jl")

# A small code to test tridiagonalization file
m = 100;
A = rand(m, m)
(P, T, Q) = tridiagonalization(A, 1 / 10^(10), m+1)
print(norm(P*T*Q' - A))
