#function usymqr()

#    return 0;

#end

include("tridiagonalization.jl")

# A small code to test tridiagonalization file
m = 10;
A = rand(m, m)
(P, T, Q, beta) = tridiagonalization(A, 1 / 10^(10), m+1)
last_beta = zeros(1,size(T,2))
last_beta[end] = beta;
#augment T to get S
S = [T ; last_beta];

#now do least squares with QR fact on S
