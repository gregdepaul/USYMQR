include("tridiagonalization.jl")

# A small code to test tridiagonalization file

A = rand(10, 10)
(P, T, Q) = tridiagonalization(A)
print(norm(P'*A*Q - T))
