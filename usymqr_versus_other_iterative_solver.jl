using IterativeSolvers
include("usymqr.jl")
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

x1, hist1 = usymqr(A, b, maxiter = 10n, tol = tol, log = true)
x2, hist2 = lsqr(A, b, atol=tol, btol=tol, conlim=1e10, maxiter=10n, log = true)
x3, hist3 = bicgstabl(A, b, 2, max_mv_products = 10n, log = true)

norm(A*x2 - b)

err1 = hist1[:resnorm];
err2 = hist2[:resnorm];
err3 = hist3[:resnorm];

using Plots
gr(reuse=true)
plot(err1, yscale=:log10)
plot!(err2, yscale=:log10)
plot!(err3, yscale=:log10)
gui()