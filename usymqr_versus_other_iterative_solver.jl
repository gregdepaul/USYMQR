using IterativeSolvers
include("usymqr.jl")
include("MarketMatrix.jl")
using Base.Test
using LinearMaps

function unsymmetric_problem(T, n)
    A = rand(T, n, n)
    x = ones(T, n)
    b = A * x
    A, x, b
end

srand(123)
n = 100

T = Float32
A, x, b = unsymmetric_problem(T, n)
#tol = sqrt(eps(real(T)))
tol = 1.0e-6


#A= Array(MatrixMarket.mmread("hydr1/hydr1.mtx"));
A= Array(MatrixMarket.mmread("g7jac010/g7jac010.mtx"));
r,c = size(A);
n = max(r,c);
x = ones(T, n);
b = A*x;

x1, hist1 = usymqr(A, b, maxiter = 1000, tol = tol, log = true)
x2, hist2 = lsqr(A, b, atol=tol, btol=tol, conlim=1e10, maxiter=1000, log = true)
x3, hist3 = bicgstabl(A, b, 2, max_mv_products = 1000, log = true)

norm(A*x2 - b)

err1 = hist1[:resnorm];
err2 = hist2[:resnorm];
err3 = hist3[:resnorm];

#labels = ["USYMQR", "LSQR", "BICGSTABL"];
labels = ["USYMQR", "LSQR"];
using Plots
gr(reuse=true)
plot(err1, yscale=:log10, title="Hollinger/g7jac010 Residual Plot", xlab="Iteration",ylab="Residual", labels = "USYMQR")
plot!(err2, yscale=:log10, labels = "LSQR")
plot!(err3, yscale=:log10, labels = "BICGSTABL")
savefig("g7jac010.png")
gui()
