#function usymqr()

#    return 0;

#end
using IterativeSolvers
include("tridiagonalization.jl")

# A small code to test tridiagonalization file
m = 10;
A = rand(m, m)
(P, T, Q, beta) = tridiagonalization(A, 1 / 10^(10), m+1)
last_beta = zeros(1,size(T,2))
last_beta[end] = beta;
#augment T to get S
S = Array([T ; last_beta]);

QR = qrfact(S);
Q = QR[:Q];
R = QR[:R];
#x = minres(S,last_beta);

#now do least squares with QR fact on S
M = Q*inv(R);
b = rand(1, size(A,1));
x0 = rand(1, size(A,1));
b1 = norm(b'-A*x0');
beta1 = zeros(1,size(S,1))';
beta1[end] = r0;
type()
x = lsqr(S,beta1);
