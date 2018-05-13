include("tridiagonalization.jl")

# A small code to test tridiagonalization file

A = rand(40, 40)
(P, T, Q) = tridiagonalization(A, 1 / 10^(4))
print(norm(P'*A*Q - T))
