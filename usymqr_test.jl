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
n = 10

T = Float32
# A, x, b = unsymmetric_problem(T, n)
# tol = sqrt(eps(real(T)))
# #tol = 1e-6;
# #x0 = rand(T, n)
# #A = Float32[[ 0.4173    0.9448    0.3377] ; [0.0497    0.4909    0.9001]; [0.9027    0.4893    0.3692]];
# #b= Float32[ 0.1112;0.7803;0.3897];
# x0 = zeros(T,n)
# x1, hist1 = usymqr(A, b, maxiter = 10n, tol = tol, log = true)
# x2, hist2 = usymqr!(x0, A, b, maxiter = 10n, tol = tol, verbose = true, log = true)
# norm(b - A * x1) / norm(b) ≤ tol
# hist1.isconverged
# isa(hist1, ConvergenceHistory)
# x2 == x0
# norm(b - A * x2) / norm(b) ≤ tol

A, x, b = unsymmetric_problem(T, n)
tol = sqrt(eps(real(T)))
x0 = rand(T, n)

x1, hist1 = usymqr(A, b, maxiter = 10n, tol = tol, log = true)
#x2, hist2 = usymqr!(x0, A, b, maxiter = 10n, tol = tol, log = true, verbose = true, log = true)
isa(hist1, ConvergenceHistory)
norm(b - A * x1) / norm(b) ≤ tol
hist1.isconverged

# @testset "Unsymmetric Matrix{$T}" for T in (Float32, Float64, Complex64, Complex128)
#     A, x, b = unsymmetric_problem(T, n)
#     tol = sqrt(eps(real(T)))
#     x0 = rand(T, n)
#
#     x1, hist1 = usymqr(A, b, maxiter = 10n, tol = tol, log = true)
#     #x2, hist2 = usymqr!(x0, A, b, maxiter = 10n, tol = tol, log = true, verbose = true, log = true)
#
#     @test isa(hist1, ConvergenceHistory)
#     @test norm(b - A * x1) / norm(b) ≤ tol
#     @test hist1.isconverged
#     #@test norm(b - A * x2) / norm(b) ≤ tol
# end

# @testset "SparseMatrixCSC{$T}" for T in (Float32, Float64, Complex64, Complex128)
#     A = let
#         B = sprand(n, n, 2 / n)
#     end
#
#     x = ones(T, n)
#     b = A * x
#     tol = sqrt(eps(real(T)))
#
#     x_approx, hist = minres(A, b, maxiter = 10n, tol = tol, log = true)
#
#     @test norm(b - A * x_approx) / norm(b) ≤ tol
#     @test hist.isconverged
# end
#
# end
